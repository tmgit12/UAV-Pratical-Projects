#!/usr/bin/env python3
"""
threat_detection_node.py
========================
Aggregates ArUco marker detections from all drones (via ros2_aruco), classifies
them as threats, landing pads, or decoys, transforms their positions to the
world frame using TF2, and publishes deduplicated maps for the swarm manager
and for scoring.

In the ICUAS '26 world (DICT_5X5_250 markers):
  • Landing pads  – markers placed ON top of buildings; landing here freezes
                    battery drain.  IDs configured via 'landing_pad_ids' param.
  • Threats       – mission targets to localize; scored by reported world pos.
                    All IDs NOT in landing_pad_ids or decoy_ids.
  • Decoys        – look like threats but incur a penalty if reported.
                    IDs configured via 'decoy_ids' param.

Subscriptions  (for each drone 1..N):
  /cf_x/aruco_markers   ros2_aruco_interfaces/msg/ArucoMarkers
  /cf_x/pose            geometry_msgs/PoseStamped   (fallback if TF unavailable)

Publications:
  /threat_detection/threats       geometry_msgs/PoseArray  (world frame)
  /threat_detection/landing_pads  geometry_msgs/PoseArray  (world frame)
  /threat_detection/all_markers   visualization_msgs/MarkerArray  (RViz)
  /threat_detection/report        std_msgs/String  (JSON — full log)

Parameters:
  num_drones        int    5
  landing_pad_ids   str    '1,2,3,4,5'   comma-separated ArUco IDs
  decoy_ids         str    ''             comma-separated ArUco IDs
  dedup_radius      float  1.0   m — same marker if detections within this
  confidence_count  int    2     — require N detections before reporting
  world_frame       str    'world'
  report_rate       float  1.0   Hz
"""

import json
import math
from collections import defaultdict
from functools import partial
from typing import Optional

import numpy as np

import rclpy
from rclpy.node import Node
from geometry_msgs.msg import PoseArray, Pose, PoseStamped
from std_msgs.msg import String
from visualization_msgs.msg import Marker, MarkerArray
from builtin_interfaces.msg import Duration

# TF2 for camera→world frame transform
from tf2_ros import Buffer, TransformListener, LookupException
from tf2_ros import ConnectivityException, ExtrapolationException

try:
    from ros2_aruco_interfaces.msg import ArucoMarkers
    ARUCO_AVAILABLE = True
except ImportError:
    ARUCO_AVAILABLE = False


# ─────────────────────────────────────────────────────────────────────────────
# Data classes
# ─────────────────────────────────────────────────────────────────────────────

class MarkerObservation:
    """Accumulates detections of a single physical marker."""

    def __init__(self, marker_id: int, world_pos: np.ndarray,
                 detected_by: int):
        self.marker_id    = marker_id
        self.positions    = [world_pos.copy()]   # running list for averaging
        self.detected_by  = {detected_by}
        self.count        = 1

    @property
    def mean_pos(self) -> np.ndarray:
        return np.mean(self.positions, axis=0)

    def update(self, world_pos: np.ndarray, detected_by: int):
        self.positions.append(world_pos.copy())
        self.detected_by.add(detected_by)
        self.count += 1


# ─────────────────────────────────────────────────────────────────────────────
# Node
# ─────────────────────────────────────────────────────────────────────────────

class ThreatDetectionNode(Node):

    def __init__(self):
        super().__init__('threat_detection_node')

        # ── Parameters ───────────────────────────────────────────
        self.declare_parameter('num_drones',       5)
        self.declare_parameter('landing_pad_ids',  '1,2,3,4,5')
        self.declare_parameter('decoy_ids',        '')
        self.declare_parameter('dedup_radius',     1.0)
        self.declare_parameter('confidence_count', 2)
        self.declare_parameter('world_frame',      'world')
        self.declare_parameter('report_rate',      1.0)

        p = self.get_parameter
        self.num_drones   = p('num_drones').value
        self.dedup_radius = p('dedup_radius').value
        self.conf_count   = p('confidence_count').value
        self.world_frame  = p('world_frame').value

        def parse_ids(param_name: str) -> set[int]:
            raw = p(param_name).value.strip()
            if not raw:
                return set()
            try:
                return {int(x.strip()) for x in raw.split(',') if x.strip()}
            except ValueError:
                self.get_logger().warn(
                    f"Could not parse '{param_name}': '{raw}'")
                return set()

        self.landing_pad_ids = parse_ids('landing_pad_ids')
        self.decoy_ids       = parse_ids('decoy_ids')

        # ── TF2 ──────────────────────────────────────────────────
        self.tf_buffer   = Buffer()
        self.tf_listener = TransformListener(self.tf_buffer, self)

        # ── Drone pose cache (fallback when TF is unavailable) ───
        self.drone_poses: dict[int, PoseStamped] = {}

        # ── Marker storage ───────────────────────────────────────
        # key: ArUco marker_id → MarkerObservation
        self.observations: dict[int, MarkerObservation] = {}

        # ── Subscriptions ────────────────────────────────────────
        for i in range(1, self.num_drones + 1):
            if ARUCO_AVAILABLE:
                self.create_subscription(
                    ArucoMarkers,
                    f'/cf_{i}/aruco_markers',
                    partial(self._aruco_cb, i),
                    10
                )
            self.create_subscription(
                PoseStamped,
                f'/cf_{i}/pose',
                partial(self._drone_pose_cb, i),
                10
            )

        # ── Publishers ───────────────────────────────────────────
        self.threats_pub  = self.create_publisher(
            PoseArray,   '/threat_detection/threats',      10)
        self.pads_pub     = self.create_publisher(
            PoseArray,   '/threat_detection/landing_pads', 10)
        self.viz_pub      = self.create_publisher(
            MarkerArray, '/threat_detection/all_markers',  10)
        self.report_pub   = self.create_publisher(
            String,      '/threat_detection/report',       10)

        rate = p('report_rate').value
        self.create_timer(1.0 / rate, self._publish_tick)

        if not ARUCO_AVAILABLE:
            self.get_logger().warn(
                "ros2_aruco_interfaces not found — "
                "ArUco detection disabled.  Only pose cache active."
            )

        self.get_logger().info(
            f"ThreatDetectionNode ready  |  "
            f"landing_pad_ids={self.landing_pad_ids}  "
            f"decoy_ids={self.decoy_ids}  "
            f"dedup_radius={self.dedup_radius} m  "
            f"confidence_count={self.conf_count}"
        )

    # ── Callbacks ──────────────────────────────────────────────────────────

    def _drone_pose_cb(self, drone_id: int, msg: PoseStamped):
        self.drone_poses[drone_id] = msg

    def _aruco_cb(self, drone_id: int, msg):
        """
        Receive ArucoMarkers from drone 'drone_id'.
        The poses in msg.poses are in the camera frame (cf_x/camera).
        Transform each to world frame via TF2, then update the observation map.
        """
        if not msg.marker_ids:
            return

        camera_frame = f'cf_{drone_id}/camera'

        for aruco_id, cam_pose in zip(msg.marker_ids, msg.poses):
            aruco_id = int(aruco_id)

            # ── Transform camera frame → world frame ─────────────
            world_pos = self._transform_to_world(
                cam_pose, camera_frame, msg.header.stamp)

            if world_pos is None:
                # TF unavailable — use drone pose as rough estimate
                world_pos = self._estimate_world_pos_from_drone(
                    drone_id, cam_pose)

            if world_pos is None:
                continue

            self.get_logger().debug(
                f"[cf_{drone_id}] Detected ArUco id={aruco_id} "
                f"world=({world_pos[0]:.2f},{world_pos[1]:.2f},{world_pos[2]:.2f})"
            )

            self._update_observations(aruco_id, world_pos, drone_id)

    # ── Coordinate transforms ──────────────────────────────────────────────

    def _transform_to_world(
        self,
        cam_pose,          # geometry_msgs/Pose  (in camera frame)
        camera_frame: str,
        stamp,
    ) -> Optional[np.ndarray]:
        """
        Use TF2 to transform a pose from the camera frame to the world frame.
        Returns np.array([x, y, z]) or None on failure.
        """
        try:
            t = self.tf_buffer.lookup_transform(
                self.world_frame, camera_frame,
                rclpy.time.Time(), timeout=rclpy.duration.Duration(seconds=0.1)
            )
        except (LookupException, ConnectivityException,
                ExtrapolationException) as e:
            self.get_logger().debug(f"TF lookup failed: {e}", throttle_duration_sec=5.0)
            return None

        # Apply transform manually (rotation + translation)
        # TF translation
        tx = t.transform.translation.x
        ty = t.transform.translation.y
        tz = t.transform.translation.z

        # Rotation quaternion
        qx = t.transform.rotation.x
        qy = t.transform.rotation.y
        qz = t.transform.rotation.z
        qw = t.transform.rotation.w

        # Point in camera frame
        px = cam_pose.position.x
        py = cam_pose.position.y
        pz = cam_pose.position.z

        # Rotate point: p_world = R * p_cam + t
        # Using quaternion rotation formula: v' = q * v * q^{-1}
        vx, vy, vz = _quat_rotate(px, py, pz, qx, qy, qz, qw)

        return np.array([vx + tx, vy + ty, vz + tz])

    def _estimate_world_pos_from_drone(
        self,
        drone_id: int,
        cam_pose,
    ) -> Optional[np.ndarray]:
        """
        Fallback when TF is unavailable.
        Approximates world position using drone pose + camera offset.
        Assumes camera looks forward along drone's +x body axis.
        """
        if drone_id not in self.drone_poses:
            return None

        dp = self.drone_poses[drone_id]
        dx = dp.pose.position.x
        dy = dp.pose.position.y
        dz = dp.pose.position.z

        q = dp.pose.orientation
        drone_yaw = math.atan2(
            2.0 * (q.w * q.z + q.x * q.y),
            1.0 - 2.0 * (q.y * q.y + q.z * q.z)
        )

        # Detected marker position relative to camera ≈ relative to drone
        rel_x = cam_pose.position.z   # camera +z is forward = body +x
        rel_y = -cam_pose.position.x  # camera +x is right = body +y
        rel_z = -cam_pose.position.y  # camera -y is up = body +z

        # Rotate by drone yaw
        wx = dx + rel_x * math.cos(drone_yaw) - rel_y * math.sin(drone_yaw)
        wy = dy + rel_x * math.sin(drone_yaw) + rel_y * math.cos(drone_yaw)
        wz = dz + rel_z

        return np.array([wx, wy, wz])

    # ── Observation map ────────────────────────────────────────────────────

    def _update_observations(
        self,
        marker_id:   int,
        world_pos:   np.ndarray,
        detected_by: int,
    ):
        """
        Update the observation map with a new detection.
        Deduplication: if a known observation of the same marker_id is within
        dedup_radius, merge; otherwise create a new entry.
        (Multiple physical markers could share the same ID in theory, but
        the competition uses unique IDs per marker.)
        """
        if marker_id in self.observations:
            obs = self.observations[marker_id]
            if np.linalg.norm(world_pos - obs.mean_pos) < self.dedup_radius:
                obs.update(world_pos, detected_by)
                return

        # New observation
        self.observations[marker_id] = MarkerObservation(
            marker_id, world_pos, detected_by)
        self.get_logger().info(
            f"New marker observed: id={marker_id}  "
            f"pos=({world_pos[0]:.2f},{world_pos[1]:.2f},{world_pos[2]:.2f})  "
            f"class={self._classify(marker_id)}"
        )

    # ── Classification ────────────────────────────────────────────────────

    def _classify(self, marker_id: int) -> str:
        """
        Classify an ArUco marker ID.
        Returns 'landing_pad', 'decoy', or 'threat'.
        """
        if marker_id in self.landing_pad_ids:
            return 'landing_pad'
        if marker_id in self.decoy_ids:
            return 'decoy'
        return 'threat'

    # ── Publication tick ──────────────────────────────────────────────────

    def _publish_tick(self):
        now_msg = self.get_clock().now().to_msg()

        threats_poses   = []
        pad_poses       = []
        marker_array    = MarkerArray()
        report_entries  = []

        for marker_id, obs in self.observations.items():
            # Only report markers with sufficient confidence
            if obs.count < self.conf_count:
                continue

            cls = self._classify(marker_id)
            pos = obs.mean_pos
            report_entries.append({
                'id':      marker_id,
                'class':   cls,
                'pos':     pos.tolist(),
                'count':   obs.count,
                'drones':  sorted(obs.detected_by),
            })

            # Build Pose
            p = Pose()
            p.position.x = float(pos[0])
            p.position.y = float(pos[1])
            p.position.z = float(pos[2])
            p.orientation.w = 1.0

            if cls == 'landing_pad':
                pad_poses.append(p)
            elif cls == 'threat':
                threats_poses.append(p)
            # decoys are logged but not published

            # RViz marker
            m = Marker()
            m.header.frame_id = self.world_frame
            m.header.stamp    = now_msg
            m.ns   = cls
            m.id   = marker_id
            m.type = Marker.CUBE
            m.action = Marker.ADD
            m.pose     = p
            m.scale.x  = m.scale.y = m.scale.z = 0.4
            m.color.a  = 0.85
            # Colour by classification
            if cls == 'threat':
                m.color.r, m.color.g, m.color.b = 1.0, 0.1, 0.1
            elif cls == 'landing_pad':
                m.color.r, m.color.g, m.color.b = 0.1, 0.8, 0.1
            else:  # decoy
                m.color.r, m.color.g, m.color.b = 1.0, 0.9, 0.0
            m.lifetime = Duration(sec=5)
            marker_array.markers.append(m)

            # Text label
            lbl = Marker()
            lbl.header  = m.header
            lbl.ns      = f'{cls}_label'
            lbl.id      = marker_id + 10000
            lbl.type    = Marker.TEXT_VIEW_FACING
            lbl.action  = Marker.ADD
            lbl.pose    = p
            lbl.pose.position.z = float(pos[2]) + 0.5
            lbl.scale.z = 0.3
            lbl.color   = m.color
            lbl.color.a = 1.0
            lbl.text    = f'{cls[:3].upper()}:{marker_id}'
            lbl.lifetime = Duration(sec=5)
            marker_array.markers.append(lbl)

        # Publish PoseArrays
        header_kwargs = dict(frame_id=self.world_frame, stamp=now_msg)

        def make_pa(poses):
            pa = PoseArray()
            pa.header.frame_id = self.world_frame
            pa.header.stamp    = now_msg
            pa.poses = poses
            return pa

        self.threats_pub.publish(make_pa(threats_poses))
        self.pads_pub.publish(make_pa(pad_poses))
        self.viz_pub.publish(marker_array)

        # JSON report
        report = {
            'threats':      [e for e in report_entries if e['class'] == 'threat'],
            'landing_pads': [e for e in report_entries if e['class'] == 'landing_pad'],
            'decoys':       [e for e in report_entries if e['class'] == 'decoy'],
        }
        self.report_pub.publish(String(data=json.dumps(report)))

        if report_entries:
            self.get_logger().info(
                f"Markers: {len(threats_poses)} threat(s)  "
                f"{len(pad_poses)} landing pad(s)  "
                f"{sum(1 for e in report_entries if e['class']=='decoy')} decoy(s)",
                throttle_duration_sec=5.0
            )


# ─────────────────────────────────────────────────────────────────────────────
# Pure-Python quaternion rotation helper
# ─────────────────────────────────────────────────────────────────────────────

def _quat_rotate(vx: float, vy: float, vz: float,
                 qx: float, qy: float, qz: float, qw: float
                 ) -> tuple[float, float, float]:
    """
    Rotate vector (vx, vy, vz) by unit quaternion (qx, qy, qz, qw).
    Uses the formula: v' = q * v * q^{-1}  where v is treated as a pure quaternion.
    """
    # Compute q * v
    wx = qw * vx + qy * vz - qz * vy
    wy = qw * vy + qz * vx - qx * vz
    wz = qw * vz + qx * vy - qy * vx
    ww = -qx * vx - qy * vy - qz * vz

    # Compute (q * v) * q^{-1}   (q^{-1} = conjugate for unit quaternion)
    rx = wx * qw + ww * (-qx) + wy * (-qz) - wz * (-qy)
    ry = wy * qw + ww * (-qy) + wz * (-qx) - wx * (-qz)
    rz = wz * qw + ww * (-qz) + wx * (-qy) - wy * (-qx)

    return rx, ry, rz


# ─────────────────────────────────────────────────────────────────────────────

def main(args=None):
    rclpy.init(args=args)
    node = ThreatDetectionNode()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == '__main__':
    main()
