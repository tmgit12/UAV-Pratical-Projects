#!/usr/bin/env python3
"""
swarm_los_network.py

Controls N Crazyflie drones to maintain a dynamic Line-of-Sight (LoS)
daisy-chain between the world origin (0,0) and a moving rover (AGV).

Drones are distributed evenly along the straight line from Origin → AGV:

    P_target[i] = P_origin + (i / (N+1)) * (P_agv - P_origin)
    for i in 1 .. N  (endpoints excluded, so drones never sit on top of
                       the origin or the AGV itself)

Key features:
  ┌──────────────────────────────────────────────────────────────┐
  │  ● Per-drone state machines (non-blocking, all run in sync)  │
  │  ● Hungarian algorithm for optimal drone↔target assignment   │
  │  ● Artificial Potential Fields (APF) for collision avoidance │
  │  ● Battery monitoring with graceful landing + network resize │
  │  ● Actual yaw from quaternion (no integration drift)         │
  │  ● cmd_hover for stable altitude independent of tilt         │
  └──────────────────────────────────────────────────────────────┘

Usage (inside the container, empty world, 2 drones):
    python3 swarm_los_network.py

With parameters:
    python3 swarm_los_network.py --ros-args \
        -p num_drones:=2 \
        -p follow_height:=1.5 \
        -p safe_distance:=1.2 \
        -p Kp_attract:=0.4 \
        -p Kr_repulse:=0.3 \
        -p battery_threshold:=20.0

ROS 2 Parameters:
    num_drones        [int,   default 2]     Number of drones (cf_1 … cf_N)
    follow_height     [float, default 1.5]   Hover altitude in metres (z_distance)
    takeoff_height    [float, default 1.0]   Initial takeoff target in metres
    takeoff_duration  [float, default 3.0]   Seconds given to the takeoff service
    height_threshold  [float, default 0.85]  Fraction of takeoff_height needed
                                             before switching to FOLLOWING
    control_rate      [float, default 20.0]  Hz for the main control loop
    max_speed         [float, default 0.5]   m/s cap on total horizontal speed
    Kp_attract        [float, default 0.4]   P-gain: distance to target → speed
    Kr_repulse        [float, default 0.3]   Repulsive gain for APF
    safe_distance     [float, default 1.2]   Min drone–drone distance (m)
                                             below which repulsion activates
    battery_threshold [float, default 20.0]  Battery % at which drone lands
    takeoff_stagger   [float, default 1.0]   Seconds between individual takeoffs
    crash_z_threshold [float, default 0.08]  z below which drone = crashed (m)
"""

import math
from enum import Enum, auto
from functools import partial

import numpy as np
from scipy.optimize import linear_sum_assignment

import rclpy
from rclpy.node import Node
from geometry_msgs.msg import Point, PoseStamped
from sensor_msgs.msg import BatteryState
from crazyflie_interfaces.msg import Hover
from crazyflie_interfaces.srv import Takeoff, Land
from builtin_interfaces.msg import Duration


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def quat_to_yaw(qx: float, qy: float, qz: float, qw: float) -> float:
    """Extract yaw (Z-rotation) from a unit quaternion — no integration drift."""
    siny_cosp = 2.0 * (qw * qz + qx * qy)
    cosy_cosp = 1.0 - 2.0 * (qy * qy + qz * qz)
    return math.atan2(siny_cosp, cosy_cosp)


def world_to_body(vx_w: float, vy_w: float, yaw: float):
    """
    Rotate a world-frame 2-D velocity vector into the drone body frame.

    Body frame is rotated by 'yaw' relative to world:
        vx_body =  vx_world · cos(yaw) + vy_world · sin(yaw)
        vy_body = -vx_world · sin(yaw) + vy_world · cos(yaw)
    """
    cy, sy = math.cos(yaw), math.sin(yaw)
    return vx_w * cy + vy_w * sy, -vx_w * sy + vy_w * cy


def to_duration(seconds: float) -> Duration:
    sec = int(seconds)
    nsec = int((seconds - sec) * 1e9)
    return Duration(sec=sec, nanosec=nsec)


# ─────────────────────────────────────────────────────────────────────────────
# Per-drone state
# ─────────────────────────────────────────────────────────────────────────────

class DroneState(Enum):
    INIT      = auto()   # waiting for services + first pose data
    TAKING_OFF = auto()  # takeoff service sent; waiting for climb
    CLIMBING  = auto()   # monitoring z until height confirmed
    FOLLOWING = auto()   # active APF control toward assigned target
    LANDING   = auto()   # battery low; land service sent
    LANDED    = auto()   # on ground, excluded from network
    CRASHED   = auto()   # unexpected z=0 while airborne


class DroneAgent:
    """
    All mutable state for one physical drone.

    The main node keeps a list of these and operates on them every tick,
    keeping the control loop fully non-blocking (no time.sleep, no spin_until).
    """

    def __init__(self, drone_id: int):
        self.id             = drone_id
        self.prefix         = f"cf_{drone_id}"

        # ── Measurements ──────────────────────────────────────────
        self.pos            = None    # np.array [x, y, z]  (world frame)
        self.yaw            = 0.0    # radians, from quaternion
        self.battery_pct    = 100.0  # 0-100 %

        # ── Control ───────────────────────────────────────────────
        self.assigned_target = None  # np.array [x, y]  (world frame, 2-D)
        self.state           = DroneState.INIT
        self.takeoff_time    = None  # rclpy.Time when takeoff service was sent

        # ── Service handles (set by the node after creation) ──────
        self.takeoff_client = None
        self.land_client    = None

        # ── Publisher handle (set by the node) ────────────────────
        self.cmd_hover_pub  = None

    # ── Convenience properties ────────────────────────────────────

    @property
    def is_active(self) -> bool:
        """True when the drone is airborne and contributing to the network."""
        return self.state in (DroneState.CLIMBING, DroneState.FOLLOWING)

    @property
    def is_controllable(self) -> bool:
        """True when we should send cmd_hover commands to this drone."""
        return self.state == DroneState.FOLLOWING

    def publish_hover(self, vx: float, vy: float, yaw_rate: float,
                      z_distance: float):
        msg            = Hover()
        msg.vx         = float(vx)
        msg.vy         = float(vy)
        msg.yaw_rate   = float(yaw_rate)
        msg.z_distance = float(z_distance)
        self.cmd_hover_pub.publish(msg)

    def hover_in_place(self, z_distance: float):
        """Stop horizontal motion; maintain altitude."""
        self.publish_hover(0.0, 0.0, 0.0, z_distance)


# ─────────────────────────────────────────────────────────────────────────────
# Main node
# ─────────────────────────────────────────────────────────────────────────────

class SwarmLOSNetwork(Node):

    ORIGIN = np.array([0.0, 0.0])   # fixed anchor of the LoS chain

    def __init__(self):
        super().__init__('swarm_los_network')

        # ── Parameters ───────────────────────────────────────────
        self.declare_parameter('num_drones',        2)
        self.declare_parameter('follow_height',     1.5)
        self.declare_parameter('takeoff_height',    1.0)
        self.declare_parameter('takeoff_duration',  3.0)
        self.declare_parameter('height_threshold',  0.85)
        self.declare_parameter('control_rate',      20.0)
        self.declare_parameter('max_speed',         0.5)
        self.declare_parameter('Kp_attract',        0.4)
        self.declare_parameter('Kr_repulse',        0.3)
        self.declare_parameter('safe_distance',     1.2)
        self.declare_parameter('battery_threshold', 20.0)
        self.declare_parameter('takeoff_stagger',   1.0)
        self.declare_parameter('crash_z_threshold', 0.08)

        p = self.get_parameter
        self.num_drones        = p('num_drones').value
        self.follow_height     = p('follow_height').value
        self.takeoff_height    = p('takeoff_height').value
        self.takeoff_duration  = p('takeoff_duration').value
        self.height_threshold  = p('height_threshold').value
        self.control_rate      = p('control_rate').value
        self.max_speed         = p('max_speed').value
        self.Kp_attract        = p('Kp_attract').value
        self.Kr_repulse        = p('Kr_repulse').value
        self.safe_distance     = p('safe_distance').value
        self.battery_threshold = p('battery_threshold').value
        self.takeoff_stagger   = p('takeoff_stagger').value
        self.crash_z_threshold = p('crash_z_threshold').value

        self.dt = 1.0 / self.control_rate

        # ── AGV state ─────────────────────────────────────────────
        self.agv_pos = None   # np.array [x, y]

        # ── Create per-drone agents ───────────────────────────────
        self.drones: list[DroneAgent] = [
            DroneAgent(i) for i in range(1, self.num_drones + 1)
        ]

        # ── Subscribers, publishers, service clients ──────────────
        self.create_subscription(
            Point, '/AGV/pose', self._agv_cb, 10)

        for drone in self.drones:
            self._setup_drone_io(drone)

        # ── Main control timer ────────────────────────────────────
        self.create_timer(self.dt, self._tick)

        self.get_logger().info(
            f"SwarmLOSNetwork ready  |  num_drones={self.num_drones}  |  "
            f"follow_height={self.follow_height} m  |  "
            f"safe_distance={self.safe_distance} m  |  "
            f"Kp={self.Kp_attract}  Kr={self.Kr_repulse}"
        )
        self.start_time = self.get_clock().now()

    # ─────────────────────────────────────────────────────────────
    # I/O setup
    # ─────────────────────────────────────────────────────────────

    def _setup_drone_io(self, drone: DroneAgent):
        """Create all ROS 2 I/O for one drone."""
        pfx = drone.prefix

        # Pose
        self.create_subscription(
            PoseStamped, f'/{pfx}/pose',
            partial(self._pose_cb, drone.id), 10)

        # Battery
        self.create_subscription(
            BatteryState, f'/{pfx}/battery_status',
            partial(self._battery_cb, drone.id), 10)

        # cmd_hover publisher
        drone.cmd_hover_pub = self.create_publisher(
            Hover, f'/{pfx}/cmd_hover', 10)

        # Service clients
        drone.takeoff_client = self.create_client(
            Takeoff, f'/{pfx}/takeoff')
        drone.land_client = self.create_client(
            Land, f'/{pfx}/land')

    # ─────────────────────────────────────────────────────────────
    # Callbacks
    # ─────────────────────────────────────────────────────────────

    def _agv_cb(self, msg: Point):
        self.agv_pos = np.array([msg.x, msg.y])

    def _pose_cb(self, drone_id: int, msg: PoseStamped):
        drone = self.drones[drone_id - 1]
        p = msg.pose.position
        drone.pos = np.array([p.x, p.y, p.z])
        q = msg.pose.orientation
        drone.yaw = quat_to_yaw(q.x, q.y, q.z, q.w)

    def _battery_cb(self, drone_id: int, msg: BatteryState):
        drone = self.drones[drone_id - 1]
        # BatteryState.percentage is in [0.0, 1.0]
        if not math.isnan(msg.percentage):
            drone.battery_pct = msg.percentage

    # ─────────────────────────────────────────────────────────────
    # Network target geometry
    # ─────────────────────────────────────────────────────────────

    def _compute_targets(self, n_active: int) -> list[np.ndarray]:
        """
        Distribute n_active waypoints evenly between Origin and AGV.

        Using i / (N+1) ensures no drone sits exactly at the origin or the
        AGV — they fill the space in between:

            P[i] = origin + (i / (N+1)) * (agv - origin),  i = 1..N
        """
        if self.agv_pos is None or n_active == 0:
            return []

        targets = []
        agv_vec = self.agv_pos - self.ORIGIN   # vector origin → AGV
        for i in range(1, n_active + 1):
            frac = i / (n_active + 1)
            targets.append(self.ORIGIN + frac * agv_vec)
        return targets

    # ─────────────────────────────────────────────────────────────
    # Hungarian assignment
    # ─────────────────────────────────────────────────────────────

    def _assign_targets(self, active_drones: list[DroneAgent],
                        targets: list[np.ndarray]):
        """
        Assign each active drone to the optimal target using the Hungarian
        algorithm (scipy linear_sum_assignment).

        Cost matrix  C[i][j] = Euclidean distance from drone i to target j.

        This prevents drones from crisscrossing when the AGV turns sharply:
        the globally optimal assignment always minimises total travel.
        """
        n = len(active_drones)
        if n == 0:
            return

        # Build n × n cost matrix
        C = np.full((n, n), 1e9)
        for i, drone in enumerate(active_drones):
            if drone.pos is None:
                continue
            for j, target in enumerate(targets):
                C[i, j] = np.linalg.norm(drone.pos[:2] - target)

        # Solve: row_ind[k] = drone index, col_ind[k] = target index
        row_ind, col_ind = linear_sum_assignment(C)
        for r, c in zip(row_ind, col_ind):
            active_drones[r].assigned_target = targets[c]

    # ─────────────────────────────────────────────────────────────
    # Artificial Potential Field (APF) velocity computation
    # ─────────────────────────────────────────────────────────────

    def _apf_velocity(self, drone: DroneAgent,
                      all_drones: list[DroneAgent]) -> tuple[float, float]:
        """
        Compute world-frame velocity for 'drone' using APF.

        ──────────────────────────────────────────────────────────
        Attractive potential:
            F_att = Kp * (target - pos_2d)          [linear P-law]
            |F_att| capped at max_speed

        Repulsive potential (for each neighbour j within safe_distance):
            d_j      = |pos_j - pos_i|
            η        = 1/d_j − 1/d_safe             (positive when d < d_safe)
            F_rep_j  = Kr * η * (1/d_j²) * unit(pos_i − pos_j)

            This is the classical APF repulsion (Khatib 1986):
            the force grows to infinity as d → 0, and vanishes at d = d_safe.

        Total:
            F_total = clip(F_att + ΣF_rep_j, 0, max_speed)
        ──────────────────────────────────────────────────────────
        """
        if drone.pos is None or drone.assigned_target is None:
            return 0.0, 0.0

        pos2 = drone.pos[:2]

        # ── Attractive force ─────────────────────────────────────
        err  = drone.assigned_target - pos2
        dist = np.linalg.norm(err)
        if dist > 0.01:
            speed = min(dist * self.Kp_attract, self.max_speed)
            F_att = (err / dist) * speed
        else:
            F_att = np.zeros(2)

        # ── Repulsive forces from neighbouring drones ─────────────
        F_rep = np.zeros(2)
        for other in all_drones:
            if other.id == drone.id or other.pos is None:
                continue
            if not other.is_active:
                continue

            diff = pos2 - other.pos[:2]     # vector pointing AWAY from other
            d    = np.linalg.norm(diff)

            if d < self.safe_distance and d > 1e-3:
                # Classical APF repulsion magnitude:
                #   η = (1/d − 1/d_safe)   — goes to ∞ as d → 0
                #   |F| = Kr * η / d²
                eta = (1.0 / d) - (1.0 / self.safe_distance)
                mag = self.Kr_repulse * eta / (d ** 2)
                F_rep += mag * (diff / d)

        # ── Combine and clip to max_speed ─────────────────────────
        F_total = F_att + F_rep
        total   = np.linalg.norm(F_total)
        if total > self.max_speed:
            F_total = F_total / total * self.max_speed

        return float(F_total[0]), float(F_total[1])

    # ─────────────────────────────────────────────────────────────
    # State machine helpers
    # ─────────────────────────────────────────────────────────────

    def _send_takeoff(self, drone: DroneAgent):
        """Fire async takeoff service call for this drone."""
        if not drone.takeoff_client.service_is_ready():
            return False
        req          = Takeoff.Request()
        req.height   = float(self.takeoff_height)
        req.duration = to_duration(self.takeoff_duration)
        future = drone.takeoff_client.call_async(req)
        future.add_done_callback(
            lambda f, d=drone: self.get_logger().info(
                f"[{d.prefix}] takeoff accepted." if not f.exception()
                else f"[{d.prefix}] takeoff FAILED: {f.exception()}"))
        drone.takeoff_time = self.get_clock().now()
        drone.state        = DroneState.CLIMBING
        self.get_logger().info(
            f"[{drone.prefix}] takeoff sent "
            f"(height={self.takeoff_height} m, "
            f"duration={self.takeoff_duration} s)")
        return True

    def _send_land(self, drone: DroneAgent):
        """Fire async land service call (battery low or crash recovery)."""
        drone.state = DroneState.LANDING
        if not drone.land_client.service_is_ready():
            self.get_logger().warn(
                f"[{drone.prefix}] land service not ready — "
                "drone will hover until service is available.")
            return
        req          = Land.Request()
        req.height   = 0.0
        req.duration = to_duration(3.0)
        future = drone.land_client.call_async(req)
        future.add_done_callback(
            lambda f, d=drone: self._landing_done(d, f))
        self.get_logger().warn(
            f"[{drone.prefix}] landing (battery={drone.battery_pct:.1f}%)")

    def _landing_done(self, drone: DroneAgent, future):
        try:
            future.result()
            self.get_logger().info(f"[{drone.prefix}] landed successfully.")
        except Exception as e:
            self.get_logger().error(f"[{drone.prefix}] land service error: {e}")
        finally:
            drone.state = DroneState.LANDED

    # ─────────────────────────────────────────────────────────────
    # Main tick
    # ─────────────────────────────────────────────────────────────

    def _tick(self):
        """
        Called at control_rate Hz.  Runs all state machines, recomputes
        targets and assignments, then publishes velocity commands.
        """

        # ── Step 1: Per-drone state transitions ──────────────────
        now          = self.get_clock().now()
        active_count = 0

        for idx, drone in enumerate(self.drones):

            # ── INIT: wait for services + first AGV pose ─────────
            if drone.state == DroneState.INIT:
                if (self.agv_pos is not None
                        and drone.takeoff_client.service_is_ready()
                        and drone.pos is not None):
                    # Stagger takeoffs: drone index * stagger_seconds
                    elapsed_s = (
                        self.get_clock().now() - self.start_time
                    ).nanoseconds * 1e-9
                    delay_done = elapsed_s > idx * self.takeoff_stagger
                    if delay_done:
                        self.get_logger().info(
                            f"[{drone.prefix}] conditions met → taking off.")
                        drone.state = DroneState.TAKING_OFF
                else:
                    missing = []
                    if self.agv_pos is None:
                        missing.append("AGV/pose")
                    if not drone.takeoff_client.service_is_ready():
                        missing.append("takeoff svc")
                    if drone.pos is None:
                        missing.append("pose")
                    self.get_logger().info(
                        f"[{drone.prefix}] INIT — waiting for: "
                        f"{', '.join(missing)}",
                        throttle_duration_sec=3.0)

            # ── TAKING_OFF: send the service call ─────────────────
            elif drone.state == DroneState.TAKING_OFF:
                self._send_takeoff(drone)

            # ── CLIMBING: wait for height confirmation ────────────
            elif drone.state == DroneState.CLIMBING:
                if drone.pos is None:
                    continue
                elapsed  = (now - drone.takeoff_time).nanoseconds * 1e-9
                min_z    = self.takeoff_height * self.height_threshold
                if (elapsed >= self.takeoff_duration
                        and drone.pos[2] >= min_z):
                    drone.state = DroneState.FOLLOWING
                    self.get_logger().info(
                        f"[{drone.prefix}] height confirmed "
                        f"(z={drone.pos[2]:.2f} m) → FOLLOWING.")

            # ── FOLLOWING: check for battery/crash ───────────────
            elif drone.state == DroneState.FOLLOWING:
                # Crash detection
                if (drone.pos is not None
                        and drone.pos[2] < self.crash_z_threshold):
                    self.get_logger().error(
                        f"[{drone.prefix}] CRASH at "
                        f"({drone.pos[0]:.2f},{drone.pos[1]:.2f},"
                        f"z={drone.pos[2]:.3f}). Removing from network.")
                    drone.state = DroneState.CRASHED
                    continue
                # Battery threshold
                if drone.battery_pct < self.battery_threshold:
                    self._send_land(drone)
                    continue

            # ── LANDING / LANDED / CRASHED: nothing to do ────────
            elif drone.state in (DroneState.LANDING,
                                 DroneState.LANDED,
                                 DroneState.CRASHED):
                pass

            if drone.is_active:
                active_count += 1

        # ── Step 2: Compute evenly-spaced target positions ───────
        #
        # Only FOLLOWING drones (not climbing) participate in the
        # assignment to avoid jitter from drones still rising.
        following = [d for d in self.drones
                     if d.state == DroneState.FOLLOWING and d.pos is not None]
        n_following = len(following)

        if n_following == 0 or self.agv_pos is None:
            return

        targets = self._compute_targets(n_following)

        # ── Step 3: Hungarian assignment ─────────────────────────
        #
        # Called every tick so the assignment adapts continuously as the
        # AGV moves.  The cost matrix is cheap for small N (2–6 drones).
        self._assign_targets(following, targets)

        # ── Step 4: APF velocity + publish cmd_hover ─────────────
        for drone in following:
            if drone.assigned_target is None:
                drone.hover_in_place(self.follow_height)
                continue

            # Compute world-frame velocity via APF
            vx_w, vy_w = self._apf_velocity(drone, self.drones)

            # Rotate to drone body frame using actual yaw from quaternion
            vx_b, vy_b = world_to_body(vx_w, vy_w, drone.yaw)

            # Publish (z_distance is constant → altitude maintained by
            # hover controller regardless of horizontal tilt)
            drone.publish_hover(vx_b, vy_b, 0.0, self.follow_height)

            # ── Logging (throttled) ───────────────────────────────
            if drone.assigned_target is not None:
                dist = np.linalg.norm(drone.assigned_target - drone.pos[:2])
                self.get_logger().info(
                    f"[{drone.prefix}] "
                    f"tgt=({drone.assigned_target[0]:.1f},"
                    f"{drone.assigned_target[1]:.1f})  "
                    f"pos=({drone.pos[0]:.1f},{drone.pos[1]:.1f},"
                    f"z={drone.pos[2]:.2f})  "
                    f"Δ={dist:.2f} m  "
                    f"bat={drone.battery_pct:.0f}%",
                    throttle_duration_sec=1.0
                )

        # ── Step 5: Network status summary ───────────────────────
        self.get_logger().info(
            f"Network: {n_following}/{self.num_drones} drones active  |  "
            f"AGV=({self.agv_pos[0]:.1f},{self.agv_pos[1]:.1f})",
            throttle_duration_sec=2.0
        )


# ─────────────────────────────────────────────────────────────────────────────
# Entry point
# ─────────────────────────────────────────────────────────────────────────────

def main(args=None):
    rclpy.init(args=args)
    node = SwarmLOSNetwork()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == '__main__':
    main()
