#!/usr/bin/env python3
"""
follow_agv.py

Makes a single Crazyflie take off and continuously hover directly above
the AGV, tracking its position in real time.

Control architecture:
  - Takeoff via /cf_x/takeoff service
  - Follow via /cf_x/cmd_hover (vx, vy, yaw_rate, z_distance)
  - z_distance is held constant → altitude independent of horizontal motion

Velocity control law (horizontal P-controller in world frame):
  error  = agv_xy - cf_xy
  speed  = clip(|error| * Kp, 0, max_speed)
  v_world = (error / |error|) * speed

  Rotation to body frame uses ACTUAL drone yaw extracted from the
  pose quaternion — not an integrated estimate (which drifts).

Parameters (--ros-args -p name:=value):
  drone_id         [int,   default 1]     which cf_<id> to control
  follow_height    [float, default 1.5]   z to hold above ground (m)
  takeoff_height   [float, default 1.0]   initial takeoff target (m)
  takeoff_duration [float, default 3.0]   seconds for takeoff service
  height_threshold [float, default 0.85]  fraction of takeoff_height needed
  control_rate     [float, default 20.0]  Hz for cmd_hover stream
  max_speed        [float, default 0.5]   m/s horizontal approach speed cap
  Kp               [float, default 0.4]   P gain: distance → speed
  yaw_tracking     [bool,  default False] rotate to face AGV travel dir
                                          (keep False unless needed — adds
                                          complexity with no benefit for pure
                                          position tracking)
  crash_z_threshold[float, default 0.08]  z below which drone = crashed (m)
"""

import math
from enum import Enum, auto

import numpy as np
import rclpy
from rclpy.node import Node
from geometry_msgs.msg import Point, PoseStamped
from crazyflie_interfaces.srv import Takeoff
from crazyflie_interfaces.msg import Hover
from builtin_interfaces.msg import Duration


# ─────────────────────────────────────────────────────────────────────────────
def _quat_to_yaw(qx: float, qy: float, qz: float, qw: float) -> float:
    """Extract yaw (rotation around Z) from a unit quaternion."""
    siny_cosp = 2.0 * (qw * qz + qx * qy)
    cosy_cosp = 1.0 - 2.0 * (qy * qy + qz * qz)
    return math.atan2(siny_cosp, cosy_cosp)


def _wrap_pi(angle: float) -> float:
    """Wrap angle to [-π, π]."""
    return (angle + math.pi) % (2.0 * math.pi) - math.pi


# ─────────────────────────────────────────────────────────────────────────────
class State(Enum):
    WAITING_FOR_SERVICE = auto()
    WAITING_FOR_AGV     = auto()
    TAKING_OFF          = auto()
    CLIMBING            = auto()
    FOLLOWING           = auto()
    CRASHED             = auto()


# ─────────────────────────────────────────────────────────────────────────────
class FollowAGVNode(Node):

    def __init__(self):
        super().__init__('follow_agv')

        # ── Parameters ───────────────────────────────────────────────
        self.declare_parameter('drone_id',           1)
        self.declare_parameter('follow_height',      1.5)
        self.declare_parameter('takeoff_height',     1.0)
        self.declare_parameter('takeoff_duration',   3.0)
        self.declare_parameter('height_threshold',   0.85)
        self.declare_parameter('control_rate',       20.0)
        self.declare_parameter('max_speed',          0.5)
        self.declare_parameter('Kp',                 0.4)
        self.declare_parameter('yaw_tracking',       False)
        self.declare_parameter('crash_z_threshold',  0.08)

        self.drone_id          = self.get_parameter('drone_id').value
        self.follow_height     = self.get_parameter('follow_height').value
        self.takeoff_height    = self.get_parameter('takeoff_height').value
        self.takeoff_duration  = self.get_parameter('takeoff_duration').value
        self.height_threshold  = self.get_parameter('height_threshold').value
        self.control_rate      = self.get_parameter('control_rate').value
        self.max_speed         = self.get_parameter('max_speed').value
        self.Kp                = self.get_parameter('Kp').value
        self.yaw_tracking      = self.get_parameter('yaw_tracking').value
        self.crash_z_threshold = self.get_parameter('crash_z_threshold').value

        self.cf_prefix = f'cf_{self.drone_id}'
        self.dt        = 1.0 / self.control_rate

        # ── State ────────────────────────────────────────────────────
        self.state           = State.WAITING_FOR_SERVICE
        self.agv_pos         = None      # np.array [x, y, z]
        self.agv_pos_prev    = None
        self.cf_pos          = None      # np.array [x, y, z]
        self.cf_yaw          = 0.0       # actual yaw from pose quaternion
        self.takeoff_sent_at = None

        # ── Subscribers ──────────────────────────────────────────────
        self.create_subscription(
            Point,       '/AGV/pose',
            self._agv_cb,     10)
        self.create_subscription(
            PoseStamped, f'/{self.cf_prefix}/pose',
            self._cf_pose_cb, 10)

        # ── Publisher ────────────────────────────────────────────────
        self.cmd_hover_pub = self.create_publisher(
            Hover, f'/{self.cf_prefix}/cmd_hover', 10)

        # ── Service client ───────────────────────────────────────────
        self.takeoff_cli = self.create_client(
            Takeoff, f'/{self.cf_prefix}/takeoff')

        # ── Main loop ────────────────────────────────────────────────
        self.create_timer(self.dt, self._tick)

        self.get_logger().info(
            f"[{self.cf_prefix}] follow_agv ready  |  "
            f"follow_height={self.follow_height} m  |  "
            f"max_speed={self.max_speed} m/s  |  "
            f"Kp={self.Kp}  |  "
            f"yaw_tracking={self.yaw_tracking}  |  "
            f"rate={self.control_rate} Hz"
        )

    # ─────────────────────────────────────────────────────────────────
    # Subscribers
    # ─────────────────────────────────────────────────────────────────

    def _agv_cb(self, msg: Point):
        self.agv_pos_prev = self.agv_pos
        self.agv_pos = np.array([msg.x, msg.y, msg.z])

    def _cf_pose_cb(self, msg: PoseStamped):
        p = msg.pose.position
        self.cf_pos = np.array([p.x, p.y, p.z])
        # Extract actual yaw from quaternion — avoids integration drift
        q = msg.pose.orientation
        self.cf_yaw = _quat_to_yaw(q.x, q.y, q.z, q.w)

    # ─────────────────────────────────────────────────────────────────
    # State machine
    # ─────────────────────────────────────────────────────────────────

    def _tick(self):

        # Crash detection (FOLLOWING only — low z during CLIMBING is normal)
        if (self.state == State.FOLLOWING
                and self.cf_pos is not None
                and self.cf_pos[2] < self.crash_z_threshold):
            self.get_logger().error(
                f"[{self.cf_prefix}] CRASH DETECTED at "
                f"({self.cf_pos[0]:.2f}, {self.cf_pos[1]:.2f}, "
                f"z={self.cf_pos[2]:.3f} m). Stopping all commands.")
            self.state = State.CRASHED
            return

        if self.state == State.WAITING_FOR_SERVICE:
            if self.takeoff_cli.service_is_ready():
                self.get_logger().info(
                    f"[{self.cf_prefix}] takeoff service online. "
                    "Waiting for AGV position...")
                self.state = State.WAITING_FOR_AGV
            else:
                self.get_logger().info(
                    f"[{self.cf_prefix}] waiting for takeoff service...",
                    throttle_duration_sec=3.0)

        elif self.state == State.WAITING_FOR_AGV:
            if self.agv_pos is not None:
                self.get_logger().info(
                    f"[{self.cf_prefix}] AGV at "
                    f"{np.round(self.agv_pos, 2)}. Taking off...")
                self.state = State.TAKING_OFF
            else:
                self.get_logger().info(
                    f"[{self.cf_prefix}] waiting for /AGV/pose...",
                    throttle_duration_sec=3.0)

        elif self.state == State.TAKING_OFF:
            self._send_takeoff()
            self.takeoff_sent_at = self.get_clock().now()
            self.state = State.CLIMBING

        elif self.state == State.CLIMBING:
            elapsed  = (self.get_clock().now() -
                        self.takeoff_sent_at).nanoseconds * 1e-9
            target_z = self.takeoff_height * self.height_threshold
            height_ok = (self.cf_pos is not None
                         and self.cf_pos[2] >= target_z
                         and elapsed >= self.takeoff_duration)
            z_str = (f"{self.cf_pos[2]:.2f} m"
                     if self.cf_pos is not None else "?")
            self.get_logger().info(
                f"[{self.cf_prefix}] climbing... t={elapsed:.1f}s z={z_str}",
                throttle_duration_sec=1.0)
            if height_ok:
                self.get_logger().info(
                    f"[{self.cf_prefix}] height confirmed at z={z_str}. "
                    "Following AGV.")
                self.state = State.FOLLOWING

        elif self.state == State.FOLLOWING:
            self._stream_hover()

        elif self.state == State.CRASHED:
            pass   # terminal — restart script to retry

    # ─────────────────────────────────────────────────────────────────
    # Takeoff
    # ─────────────────────────────────────────────────────────────────

    def _send_takeoff(self):
        req          = Takeoff.Request()
        req.height   = float(self.takeoff_height)
        req.duration = _to_duration(self.takeoff_duration)
        fut = self.takeoff_cli.call_async(req)
        fut.add_done_callback(self._takeoff_cb)
        self.get_logger().info(
            f"[{self.cf_prefix}] takeoff sent "
            f"(height={self.takeoff_height} m, "
            f"duration={self.takeoff_duration} s)")

    def _takeoff_cb(self, future):
        try:
            future.result()
        except Exception as exc:
            self.get_logger().error(
                f"[{self.cf_prefix}] takeoff service failed: {exc}")

    # ─────────────────────────────────────────────────────────────────
    # Hover streaming
    # ─────────────────────────────────────────────────────────────────

    def _stream_hover(self):
        if self.agv_pos is None:
            self._publish_hover(0.0, 0.0, 0.0)
            self.get_logger().warn(
                f"[{self.cf_prefix}] AGV/pose lost — holding.",
                throttle_duration_sec=2.0)
            return
        if self.cf_pos is None:
            return

        # ── World-frame P-controller ──────────────────────────────
        err_x  = self.agv_pos[0] - self.cf_pos[0]
        err_y  = self.agv_pos[1] - self.cf_pos[1]
        dist   = math.sqrt(err_x**2 + err_y**2)

        if dist > 0.01:
            speed    = min(dist * self.Kp, self.max_speed)
            vx_world = (err_x / dist) * speed
            vy_world = (err_y / dist) * speed
        else:
            vx_world = 0.0
            vy_world = 0.0

        # ── Rotate world → body frame using ACTUAL drone yaw ──────
        #
        # We use cf_yaw extracted from the pose quaternion, not an
        # integrated estimate, so it never drifts.
        #
        # body_x =  world_x * cos(yaw) + world_y * sin(yaw)
        # body_y = -world_x * sin(yaw) + world_y * cos(yaw)
        #
        cy = math.cos(self.cf_yaw)
        sy = math.sin(self.cf_yaw)
        vx_body =  vx_world * cy + vy_world * sy
        vy_body = -vx_world * sy + vy_world * cy

        # ── Optional yaw tracking: face AGV travel direction ──────
        yaw_rate = 0.0
        if self.yaw_tracking and self.agv_pos_prev is not None:
            dx = self.agv_pos[0] - self.agv_pos_prev[0]
            dy = self.agv_pos[1] - self.agv_pos_prev[1]
            if dx * dx + dy * dy > 1e-6:
                desired_yaw = math.atan2(dy, dx)
                yaw_err     = _wrap_pi(desired_yaw - self.cf_yaw)
                yaw_rate    = float(np.clip(yaw_err * 1.0, -0.5, 0.5))

        self._publish_hover(vx_body, vy_body, yaw_rate)

        self.get_logger().info(
            f"[{self.cf_prefix}] "
            f"AGV=({self.agv_pos[0]:6.2f},{self.agv_pos[1]:6.2f})  "
            f"CF=({self.cf_pos[0]:6.2f},{self.cf_pos[1]:6.2f},"
            f"z={self.cf_pos[2]:5.2f})  "
            f"Δxy={dist:.2f} m  "
            f"yaw={math.degrees(self.cf_yaw):5.1f}°  "
            f"v_body=({vx_body:.2f},{vy_body:.2f})",
            throttle_duration_sec=1.0
        )

    def _publish_hover(self, vx: float, vy: float, yaw_rate: float):
        msg            = Hover()
        msg.vx         = float(vx)
        msg.vy         = float(vy)
        msg.yaw_rate   = float(yaw_rate)
        msg.z_distance = float(self.follow_height)
        self.cmd_hover_pub.publish(msg)


# ─────────────────────────────────────────────────────────────────────────────
def _to_duration(seconds: float) -> Duration:
    sec  = int(seconds)
    nsec = int((seconds - sec) * 1e9)
    return Duration(sec=sec, nanosec=nsec)


def main(args=None):
    rclpy.init(args=args)
    node = FollowAGVNode()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == '__main__':
    main()
