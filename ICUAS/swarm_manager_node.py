#!/usr/bin/env python3
"""
swarm_manager_node.py
=====================
Central coordinator for the ICUAS '26 swarm.  Manages up to 5 Crazyflie
drones through a per-drone state machine and assigns each one a dynamic role:

  RELAY    – Anchored at a corner waypoint in the LoS chain to the AGV.
  EXPLORER – Sweeping the map for ArUco markers when not needed for relay.

Key capabilities:
  • Non-blocking takeoff / land service calls via call_async.
  • Hungarian optimal assignment: relay drones → anchor waypoints.
  • Three-message handover protocol (REQUEST → ACCEPT → ARRIVED) ensures
    a replacement arrives before the low-battery drone leaves its post.
  • Graceful degradation: explorers return home when battery is low;
    relay drones are prioritised over exploration when relays are scarce.
  • Perching on ArUco landing pads to freeze battery drain.

Subscriptions:
  /cf_x/pose               geometry_msgs/PoseStamped
  /cf_x/battery_status     sensor_msgs/BatteryState
  /path_planning/relay_waypoints  geometry_msgs/PoseArray
  /threat_detection/landing_pads  geometry_msgs/PoseArray
  /swarm/handover          std_msgs/String  (JSON protocol bus)

Publications:
  /cf_x/cmd_hover          crazyflie_interfaces/Hover
  /swarm/handover          std_msgs/String
  /swarm/status            std_msgs/String  (JSON dashboard)

Services called:
  /cf_x/takeoff            crazyflie_interfaces/srv/Takeoff
  /cf_x/land               crazyflie_interfaces/srv/Land

Parameters:
  num_drones          int    5
  follow_height       float  1.5   m
  takeoff_height      float  1.0   m
  takeoff_duration    float  3.0   s
  height_threshold    float  0.85  fraction of takeoff_height to confirm climb
  control_rate        float  20.0  Hz
  max_speed           float  0.5   m/s
  Kp_attract          float  0.4   P-gain distance → speed
  Kr_repulse          float  0.3   APF repulsion gain
  safe_distance       float  1.2   m  inter-drone repulsion radius
  battery_warning     float  30.0  % → initiate handover
  battery_emergency   float  15.0  % → emergency land regardless
  handover_radius     float  1.5   m  replacement considered "arrived"
  takeoff_stagger     float  1.0   s  delay between successive takeoffs
  crash_z_threshold   float  0.08  m  z below which drone = crashed
  explorer_alt        float  2.5   m  explorers fly a bit higher for cameras
  search_area_half    float  20.0  m  half-side of the search square
"""

import json
import math
from enum import Enum, auto
from functools import partial
from typing import Optional

import numpy as np
from scipy.optimize import linear_sum_assignment

import rclpy
import rclpy.time
from rclpy.node import Node
from geometry_msgs.msg import Point, PoseArray, PoseStamped
from sensor_msgs.msg import BatteryState
from std_msgs.msg import String
from crazyflie_interfaces.msg import Hover
from crazyflie_interfaces.srv import Takeoff, Land
from builtin_interfaces.msg import Duration


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def quat_to_yaw(qx, qy, qz, qw) -> float:
    return math.atan2(2.0 * (qw * qz + qx * qy),
                      1.0 - 2.0 * (qy * qy + qz * qz))


def world_to_body(vx_w: float, vy_w: float, yaw: float):
    cy, sy = math.cos(yaw), math.sin(yaw)
    return vx_w * cy + vy_w * sy, -vx_w * sy + vy_w * cy


def to_duration(seconds: float) -> Duration:
    s = int(seconds)
    return Duration(sec=s, nanosec=int((seconds - s) * 1e9))


def dist2d(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a[:2] - b[:2]))


# ─────────────────────────────────────────────────────────────────────────────
# State machine
# ─────────────────────────────────────────────────────────────────────────────

class Role(Enum):
    IDLE              = auto()   # on ground, unassigned
    TAKING_OFF        = auto()   # ascending
    CLIMBING          = auto()   # waiting for height confirmation
    RELAY             = auto()   # holding anchor waypoint in LoS chain
    RELAY_MOVING      = auto()   # transitioning to a new anchor
    EXPLORER          = auto()   # sweeping for ArUco markers
    HANDOVER_OUT      = auto()   # waiting for replacement (still holding post)
    HANDOVER_IN       = auto()   # flying to take over a relay post
    RETURNING         = auto()   # low battery → flying to base
    CHARGING          = auto()   # on ground at base, recharging
    PERCHING          = auto()   # flying to ArUco landing pad
    PERCHED           = auto()   # on landing pad (battery frozen)
    EMERGENCY         = auto()   # very low battery — land NOW
    CRASHED           = auto()   # detected z < crash_threshold


# ─────────────────────────────────────────────────────────────────────────────
# Per-drone data
# ─────────────────────────────────────────────────────────────────────────────

class DroneAgent:
    """All mutable state for one Crazyflie drone."""

    # search grid step for EXPLORER lawnmower pattern (metres)
    SEARCH_STEP = 8.0

    def __init__(self, drone_id: int, search_area_half: float):
        self.id     = drone_id
        self.prefix = f'cf_{drone_id}'

        # ── Measurements ────────────────────────────────────────
        self.pos:         Optional[np.ndarray] = None   # [x, y, z]
        self.yaw:         float                = 0.0
        self.battery_pct: float                = 100.0

        # ── Role state ─────────────────────────────────────────
        self.role:          Role               = Role.IDLE
        self.target:        Optional[np.ndarray] = None  # world [x, y, z]
        self.chain_idx:     int                = -1      # relay position index

        # ── Timing ─────────────────────────────────────────────
        self.takeoff_time:  Optional[rclpy.time.Time] = None

        # ── Handover ───────────────────────────────────────────
        self.handover_partner: Optional[int] = None  # ID of other drone

        # ── Explorer search state ───────────────────────────────
        self._search_half  = search_area_half
        self._search_pts   = self._build_search_grid()
        self._search_idx   = (drone_id - 1) % max(1, len(self._search_pts))

        # ── ROS handles (set by node) ───────────────────────────
        self.cmd_hover_pub  = None
        self.takeoff_client = None
        self.land_client    = None

    # ── Properties ────────────────────────────────────────────────

    @property
    def is_airborne(self) -> bool:
        return self.role not in (Role.IDLE, Role.CHARGING, Role.PERCHED,
                                 Role.CRASHED)

    @property
    def is_relay_capable(self) -> bool:
        """Can this drone accept a RELAY assignment right now?"""
        return self.role in (Role.EXPLORER, Role.CHARGING, Role.IDLE,
                             Role.PERCHED)

    @property
    def next_search_target(self) -> np.ndarray:
        """Return next search waypoint and advance index."""
        if not self._search_pts:
            return np.zeros(3)
        pt = self._search_pts[self._search_idx % len(self._search_pts)]
        self._search_idx = (self._search_idx + 1) % len(self._search_pts)
        return pt

    def _build_search_grid(self) -> list[np.ndarray]:
        """Lawnmower pattern within [-half, half]²."""
        half = self._search_half
        step = self.SEARCH_STEP
        xs = list(np.arange(-half, half, step))
        ys = list(np.arange(-half, half, step))
        pts = []
        for i, y in enumerate(ys):
            row_x = xs if i % 2 == 0 else xs[::-1]
            for x in row_x:
                pts.append(np.array([x, y, 0.0]))  # z set at runtime
        return pts

    # ── Publishing ────────────────────────────────────────────────

    def publish_hover(self, vx: float, vy: float, yaw_rate: float,
                      z_distance: float):
        msg            = Hover()
        msg.vx         = float(vx)
        msg.vy         = float(vy)
        msg.yaw_rate   = float(yaw_rate)
        msg.z_distance = float(z_distance)
        self.cmd_hover_pub.publish(msg)

    def hover_in_place(self, z_distance: float):
        self.publish_hover(0.0, 0.0, 0.0, z_distance)


# ─────────────────────────────────────────────────────────────────────────────
# Swarm manager node
# ─────────────────────────────────────────────────────────────────────────────

class SwarmManagerNode(Node):

    BASE = np.array([0.0, 0.0, 0.0])   # base station world position

    def __init__(self):
        super().__init__('swarm_manager_node')

        # ── Parameters ───────────────────────────────────────────
        self.declare_parameter('num_drones',        5)
        self.declare_parameter('follow_height',     1.5)
        self.declare_parameter('takeoff_height',    1.0)
        self.declare_parameter('takeoff_duration',  3.0)
        self.declare_parameter('height_threshold',  0.85)
        self.declare_parameter('control_rate',      20.0)
        self.declare_parameter('max_speed',         0.4)
        self.declare_parameter('Kp_attract',        0.4)
        self.declare_parameter('Kr_repulse',        0.3)
        self.declare_parameter('safe_distance',     1.2)
        self.declare_parameter('battery_warning',   30.0)
        self.declare_parameter('battery_emergency', 15.0)
        self.declare_parameter('handover_radius',   1.5)
        self.declare_parameter('takeoff_stagger',   1.0)
        self.declare_parameter('crash_z_threshold', 0.08)
        self.declare_parameter('explorer_alt',      2.5)
        self.declare_parameter('search_area_half',  20.0)

        p = self.get_parameter
        self.num_drones        = p('num_drones').value
        self.follow_height     = p('follow_height').value
        self.takeoff_height    = p('takeoff_height').value
        self.takeoff_duration  = p('takeoff_duration').value
        self.height_threshold  = p('height_threshold').value
        self.control_rate      = p('control_rate').value
        self.max_speed         = p('max_speed').value
        self.Kp                = p('Kp_attract').value
        self.Kr                = p('Kr_repulse').value
        self.safe_dist         = p('safe_distance').value
        self.bat_warn          = p('battery_warning').value
        self.bat_emerg         = p('battery_emergency').value
        self.handover_r        = p('handover_radius').value
        self.stagger           = p('takeoff_stagger').value
        self.crash_z           = p('crash_z_threshold').value
        self.explorer_alt      = p('explorer_alt').value
        search_half            = p('search_area_half').value

        self.dt = 1.0 / self.control_rate

        # ── Drone agents ─────────────────────────────────────────
        self.drones: list[DroneAgent] = [
            DroneAgent(i, search_half)
            for i in range(1, self.num_drones + 1)
        ]

        # ── Shared world state ───────────────────────────────────
        self.relay_waypoints: list[np.ndarray] = []
        self.landing_pads:    list[np.ndarray] = []
        self.start_time = self.get_clock().now()

        # ── ROS I/O ──────────────────────────────────────────────
        for drone in self.drones:
            self._setup_io(drone)

        self.create_subscription(PoseArray, '/path_planning/relay_waypoints',
                                 self._relay_wp_cb, 10)
        self.create_subscription(PoseArray, '/threat_detection/landing_pads',
                                 self._landing_pad_cb, 10)

        self.handover_pub = self.create_publisher(String, '/swarm/handover', 10)
        self.status_pub   = self.create_publisher(String, '/swarm/status',   10)
        self.create_subscription(String, '/swarm/handover',
                                 self._handover_cb, 10)

        self.create_timer(self.dt,  self._tick)
        self.create_timer(2.0,      self._publish_status)

        self.get_logger().info(
            f"SwarmManagerNode ready  |  num_drones={self.num_drones}  |  "
            f"bat_warn={self.bat_warn}%  bat_emerg={self.bat_emerg}%"
        )

    # ── I/O setup ──────────────────────────────────────────────────────────

    def _setup_io(self, d: DroneAgent):
        pfx = d.prefix
        self.create_subscription(
            PoseStamped, f'/{pfx}/pose', partial(self._pose_cb, d.id), 10)
        self.create_subscription(
            BatteryState, f'/{pfx}/battery_status',
            partial(self._battery_cb, d.id), 10)
        d.cmd_hover_pub = self.create_publisher(Hover, f'/{pfx}/cmd_hover', 10)
        d.takeoff_client = self.create_client(Takeoff, f'/{pfx}/takeoff')
        d.land_client    = self.create_client(Land,    f'/{pfx}/land')

    # ── Callbacks ──────────────────────────────────────────────────────────

    def _pose_cb(self, drone_id: int, msg: PoseStamped):
        d = self.drones[drone_id - 1]
        p = msg.pose.position
        d.pos = np.array([p.x, p.y, p.z])
        q = msg.pose.orientation
        d.yaw = quat_to_yaw(q.x, q.y, q.z, q.w)

    def _battery_cb(self, drone_id: int, msg: BatteryState):
        d = self.drones[drone_id - 1]
        if not math.isnan(msg.percentage):
            d.battery_pct = float(msg.percentage)   # already 0-100

    def _relay_wp_cb(self, msg: PoseArray):
        """New relay anchor positions from path_planning_node."""
        self.relay_waypoints = [
            np.array([p.position.x, p.position.y, p.position.z])
            for p in msg.poses
        ]
        self._assign_relay_roles()

    def _landing_pad_cb(self, msg: PoseArray):
        self.landing_pads = [
            np.array([p.position.x, p.position.y, p.position.z])
            for p in msg.poses
        ]

    def _handover_cb(self, msg: String):
        """Process handover protocol JSON messages."""
        try:
            data = json.loads(msg.data)
        except json.JSONDecodeError:
            return

        t = data.get('type')

        if t == 'REQUEST':
            self._handle_handover_request(data)
        elif t == 'ACCEPT':
            self._handle_handover_accept(data)
        elif t == 'ARRIVED':
            self._handle_handover_arrived(data)

    # ── Role assignment (Hungarian) ────────────────────────────────────────

    def _assign_relay_roles(self):
        """
        Assign relay drones to anchor waypoints using the Hungarian algorithm
        to minimise total travel distance and avoid crisscrossing.

        Only drones already in RELAY or RELAY_MOVING roles, plus available
        idle/explorer drones, participate.  Perched/charging drones are used
        as a last resort.
        """
        n_wps = len(self.relay_waypoints)
        if n_wps == 0:
            # No relays needed — promote all airborne drones to EXPLORER
            for d in self.drones:
                if d.role == Role.RELAY:
                    d.role   = Role.EXPLORER
                    d.target = None
            return

        # Candidate drones: current relays first, then explorers, then idle
        priority_order = (
            [d for d in self.drones if d.role in (Role.RELAY, Role.RELAY_MOVING)]
            + [d for d in self.drones if d.role == Role.EXPLORER]
            + [d for d in self.drones if d.role == Role.IDLE]
        )
        candidates = priority_order[:n_wps]
        n_cand     = len(candidates)

        if n_cand == 0:
            return

        # Build cost matrix (n_cand × n_wps)
        n = max(n_cand, n_wps)
        C = np.full((n, n), 1e9)
        for i, drone in enumerate(candidates):
            if drone.pos is None:
                continue
            for j, wp in enumerate(self.relay_waypoints):
                C[i, j] = dist2d(drone.pos, wp)

        row_idx, col_idx = linear_sum_assignment(C)
        assigned = set()
        for r, c in zip(row_idx, col_idx):
            if r >= n_cand or c >= n_wps:
                continue
            drone = candidates[r]
            wp    = self.relay_waypoints[c]
            drone.target    = wp.copy()
            drone.chain_idx = c
            if drone.role not in (Role.RELAY, Role.RELAY_MOVING):
                drone.role = Role.RELAY_MOVING
            assigned.add(drone.id)

        # Drones not assigned to a relay slot become explorers
        for d in self.drones:
            if d.id not in assigned and d.role in (Role.RELAY, Role.RELAY_MOVING):
                d.role      = Role.EXPLORER
                d.chain_idx = -1
                d.target    = None

    # ── Handover protocol ─────────────────────────────────────────────────

    def _request_handover(self, dying: DroneAgent):
        """
        Stage 1: dying RELAY drone announces it needs a replacement.
        It stays at its post (HANDOVER_OUT) until replacement arrives.
        """
        dying.role = Role.HANDOVER_OUT
        msg = json.dumps({
            'type':      'REQUEST',
            'drone_id':  dying.id,
            'chain_idx': dying.chain_idx,
            'pos':       dying.pos[:3].tolist(),
            'battery':   dying.battery_pct,
        })
        self.handover_pub.publish(String(data=msg))
        self.get_logger().warn(
            f"[{dying.prefix}] Handover REQUEST  "
            f"(chain_idx={dying.chain_idx}, bat={dying.battery_pct:.0f}%)"
        )

    def _handle_handover_request(self, data: dict):
        """
        Stage 1 receiver: find the best available drone to accept the handover.
        Best = highest battery among drones that are not already RELAY.
        """
        dying_id = data['drone_id']
        if dying_id == self._own_id_guard():
            return  # we are the dying drone

        # Find best candidate
        candidates = [
            d for d in self.drones
            if d.id != dying_id
            and d.is_relay_capable
            and d.pos is not None
            and d.battery_pct > self.bat_warn + 5  # must have margin
        ]
        if not candidates:
            return

        best = max(candidates, key=lambda d: d.battery_pct)

        # Accept
        best.role          = Role.HANDOVER_IN
        best.target        = np.array(data['pos'])
        best.chain_idx     = data['chain_idx']
        best.handover_partner = dying_id

        msg = json.dumps({
            'type':       'ACCEPT',
            'from_drone': best.id,
            'to_drone':   dying_id,
            'chain_idx':  data['chain_idx'],
        })
        self.handover_pub.publish(String(data=msg))
        self.get_logger().info(
            f"[{best.prefix}] Handover ACCEPT → replacing [{dying_id}] "
            f"at chain_idx={data['chain_idx']}"
        )

    def _handle_handover_accept(self, data: dict):
        """Stage 2: dying drone records who accepted."""
        dying_id = data['to_drone']
        d = self._get_drone(dying_id)
        if d is not None and d.role == Role.HANDOVER_OUT:
            d.handover_partner = data['from_drone']

    def _handle_handover_arrived(self, data: dict):
        """
        Stage 3: replacement has arrived at the post.
        The dying drone may now leave.
        """
        dying_id = data['to_drone']
        d = self._get_drone(dying_id)
        if d is None or d.role != Role.HANDOVER_OUT:
            return

        # Dying drone: choose landing pad or base
        d.role   = Role.RETURNING
        d.target = self._nearest_landing_pad(d.pos) or self.BASE.copy()
        self.get_logger().info(
            f"[{d.prefix}] Replacement ARRIVED — departing for "
            f"({'landing pad' if self.landing_pads else 'base'})"
        )

    def _own_id_guard(self) -> int:
        """Placeholder — each agent checks its own ID in callbacks."""
        return -1

    # ── Tick ──────────────────────────────────────────────────────────────

    def _tick(self):
        now     = self.get_clock().now()
        elapsed = (now - self.start_time).nanoseconds * 1e-9

        for idx, d in enumerate(self.drones):

            # ── CRASH detection ──────────────────────────────────
            if (d.role not in (Role.IDLE, Role.CHARGING, Role.PERCHED,
                               Role.CRASHED, Role.EMERGENCY)
                    and d.pos is not None
                    and d.pos[2] < self.crash_z
                    and d.role not in (Role.TAKING_OFF, Role.CLIMBING)):
                self.get_logger().error(
                    f"[{d.prefix}] CRASH detected (z={d.pos[2]:.3f})")
                d.role = Role.CRASHED
                continue

            # ── EMERGENCY battery ─────────────────────────────────
            if (d.battery_pct < self.bat_emerg
                    and d.role not in (Role.IDLE, Role.CHARGING, Role.PERCHED,
                                       Role.CRASHED, Role.EMERGENCY)):
                self.get_logger().error(
                    f"[{d.prefix}] EMERGENCY — bat={d.battery_pct:.0f}%")
                d.role   = Role.EMERGENCY
                d.target = self.BASE.copy()
                self._send_land(d)
                continue

            # ── Battery warning for RELAY drones ──────────────────
            if (d.role == Role.RELAY
                    and d.battery_pct < self.bat_warn):
                self._request_handover(d)
                continue

            # ── Battery warning for EXPLORER drones ───────────────
            if (d.role == Role.EXPLORER
                    and d.battery_pct < self.bat_warn):
                d.role   = Role.RETURNING
                d.target = self._nearest_landing_pad(d.pos) or self.BASE.copy()
                continue

            # ── State transitions ─────────────────────────────────
            if d.role == Role.IDLE:
                if (d.pos is not None
                        and d.takeoff_client.service_is_ready()
                        and elapsed > idx * self.stagger):
                    self._send_takeoff(d, now)

            elif d.role == Role.CLIMBING:
                self._check_climb_complete(d, now)

            elif d.role == Role.RELAY:
                pass  # navigation handled below

            elif d.role == Role.RELAY_MOVING:
                if self._arrived(d, d.target):
                    d.role = Role.RELAY

            elif d.role == Role.HANDOVER_IN:
                if self._arrived(d, d.target):
                    # Announce arrival
                    msg = json.dumps({
                        'type':       'ARRIVED',
                        'from_drone': d.id,
                        'to_drone':   d.handover_partner,
                        'chain_idx':  d.chain_idx,
                    })
                    self.handover_pub.publish(String(data=msg))
                    d.role = Role.RELAY
                    self.get_logger().info(
                        f"[{d.prefix}] Now RELAY at chain_idx={d.chain_idx}")

            elif d.role == Role.EXPLORER:
                if d.target is None or self._arrived(d, d.target):
                    pt = d.next_search_target.copy()
                    pt[2] = self.explorer_alt
                    d.target = pt

            elif d.role == Role.RETURNING:
                if self._arrived_2d(d, d.target, radius=2.0):
                    self._send_land(d)
                    d.role = Role.CHARGING

            elif d.role == Role.PERCHING:
                if self._arrived_2d(d, d.target, radius=0.5):
                    self._send_land(d)
                    d.role = Role.PERCHED

            # ── Navigation commands ────────────────────────────────
            if d.role in (Role.RELAY, Role.RELAY_MOVING, Role.HANDOVER_IN,
                          Role.HANDOVER_OUT, Role.EXPLORER, Role.RETURNING,
                          Role.PERCHING, Role.EMERGENCY) and d.pos is not None:
                self._navigate(d)

    # ── Navigation: APF P-controller with body-frame rotation ─────────────

    def _navigate(self, d: DroneAgent):
        if d.target is None:
            d.hover_in_place(self.follow_height)
            return

        # Desired altitude: explorers fly higher
        z_tgt = (self.explorer_alt if d.role == Role.EXPLORER
                 else self.follow_height)

        err_x = d.target[0] - d.pos[0]
        err_y = d.target[1] - d.pos[1]
        dist  = math.sqrt(err_x ** 2 + err_y ** 2)

        # ── Attractive velocity (P-control, capped) ──────────────
        if dist > 0.01:
            speed = min(dist * self.Kp, self.max_speed)
            vx_w  = (err_x / dist) * speed
            vy_w  = (err_y / dist) * speed
        else:
            vx_w = vy_w = 0.0

        # ── Repulsive velocity (APF inter-drone) ─────────────────
        rx_w = ry_w = 0.0
        for other in self.drones:
            if other.id == d.id or other.pos is None or not other.is_airborne:
                continue
            diff_x = d.pos[0] - other.pos[0]
            diff_y = d.pos[1] - other.pos[1]
            sep    = math.sqrt(diff_x ** 2 + diff_y ** 2)
            if 1e-3 < sep < self.safe_dist:
                eta    = (1.0 / sep) - (1.0 / self.safe_dist)
                mag    = self.Kr * eta / (sep ** 2)
                rx_w  += mag * diff_x / sep
                ry_w  += mag * diff_y / sep

        vx_w += rx_w
        vy_w += ry_w

        # Clip total speed
        spd = math.sqrt(vx_w ** 2 + vy_w ** 2)
        if spd > self.max_speed:
            vx_w = vx_w / spd * self.max_speed
            vy_w = vy_w / spd * self.max_speed

        # ── Rotate to body frame using actual yaw ────────────────
        vx_b, vy_b = world_to_body(vx_w, vy_w, d.yaw)
        d.publish_hover(vx_b, vy_b, 0.0, z_tgt)

        self.get_logger().info(
            f"[{d.prefix}|{d.role.name[:4]}] "
            f"tgt=({d.target[0]:.1f},{d.target[1]:.1f}) "
            f"Δxy={dist:.2f} m  bat={d.battery_pct:.0f}%",
            throttle_duration_sec=1.0
        )

    # ── Takeoff / land ────────────────────────────────────────────────────

    def _send_takeoff(self, d: DroneAgent, now: rclpy.time.Time):
        if not d.takeoff_client.service_is_ready():
            return
        req            = Takeoff.Request()
        req.height     = float(self.takeoff_height)
        req.duration   = to_duration(self.takeoff_duration)
        future = d.takeoff_client.call_async(req)
        future.add_done_callback(lambda f, drone=d: self.get_logger().info(
            f"[{drone.prefix}] takeoff accepted" if not f.exception()
            else f"[{drone.prefix}] takeoff FAILED: {f.exception()}"))
        d.role         = Role.CLIMBING
        d.takeoff_time = now
        self.get_logger().info(f"[{d.prefix}] Takeoff sent")

    def _send_land(self, d: DroneAgent):
        if not d.land_client.service_is_ready():
            return
        req          = Land.Request()
        req.height   = 0.0
        req.duration = to_duration(3.0)
        future = d.land_client.call_async(req)
        future.add_done_callback(lambda f, drone=d: self.get_logger().info(
            f"[{drone.prefix}] land accepted" if not f.exception()
            else f"[{drone.prefix}] land FAILED: {f.exception()}"))

    def _check_climb_complete(self, d: DroneAgent, now: rclpy.time.Time):
        if d.pos is None or d.takeoff_time is None:
            return
        elapsed  = (now - d.takeoff_time).nanoseconds * 1e-9
        height_ok = (d.pos[2] >= self.takeoff_height * self.height_threshold
                     and elapsed >= self.takeoff_duration)
        if height_ok:
            self.get_logger().info(
                f"[{d.prefix}] Climb confirmed (z={d.pos[2]:.2f} m)"
            )
            # Assign role
            n_relays = len(self.relay_waypoints)
            n_existing_relays = sum(
                1 for x in self.drones
                if x.id != d.id and x.role in (Role.RELAY, Role.RELAY_MOVING))
            if n_existing_relays < n_relays:
                d.role = Role.RELAY_MOVING
                # Will get target from next _assign_relay_roles call
            else:
                d.role = Role.EXPLORER

    # ── Helpers ───────────────────────────────────────────────────────────

    def _arrived(self, d: DroneAgent, target: Optional[np.ndarray],
                 radius: float = 1.0) -> bool:
        if d.pos is None or target is None:
            return False
        return dist2d(d.pos, target) < radius

    def _arrived_2d(self, d: DroneAgent, target: Optional[np.ndarray],
                    radius: float = 1.5) -> bool:
        return self._arrived(d, target, radius)

    def _get_drone(self, drone_id: int) -> Optional[DroneAgent]:
        if 1 <= drone_id <= self.num_drones:
            return self.drones[drone_id - 1]
        return None

    def _nearest_landing_pad(self,
                              pos: Optional[np.ndarray]) -> Optional[np.ndarray]:
        if not self.landing_pads or pos is None:
            return None
        return min(self.landing_pads, key=lambda p: dist2d(pos, p)).copy()

    # ── Status publisher ──────────────────────────────────────────────────

    def _publish_status(self):
        status = {
            d.prefix: {
                'role':    d.role.name,
                'battery': round(d.battery_pct, 1),
                'chain':   d.chain_idx,
                'pos':     d.pos[:3].tolist() if d.pos is not None else None,
            }
            for d in self.drones
        }
        self.status_pub.publish(String(data=json.dumps(status)))


# ─────────────────────────────────────────────────────────────────────────────

def main(args=None):
    rclpy.init(args=args)
    node = SwarmManagerNode()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == '__main__':
    main()
