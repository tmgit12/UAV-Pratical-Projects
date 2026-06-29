#!/usr/bin/env python3
"""
path_planning_node.py
=====================
Computes an obstacle-aware path from the base station (origin) to the
moving AGV using A* on the 2D projected OccupancyGrid published by the
octomap_server.  From the raw A* path it extracts the minimum set of
"corner anchor" waypoints — positions where relay drones should park
themselves to maintain unobstructed Line-of-Sight around buildings.

Subscriptions:
  /projected_map        nav_msgs/OccupancyGrid   (from octomap_server)
  /AGV/pose             geometry_msgs/Point

Publications:
  /path_planning/relay_waypoints   geometry_msgs/PoseArray   anchor positions
  /path_planning/agv_path          nav_msgs/Path             full A* path
  /path_planning/viz               visualization_msgs/MarkerArray  (RViz)

Parameters:
  flight_height    [float, default 1.5]   z for all published waypoints (m)
  safety_margin    [float, default 0.5]   obstacle inflation radius (m)
  comm_range       [float, default 70.0]  max LoS range between nodes (m)
  replan_dist      [float, default 2.0]   replan when AGV moves this far (m)
  base_x / base_y  [float, default 0.0]  world-frame origin of base station
"""

import heapq
import math
import json
from typing import Optional

import numpy as np
import rclpy
from rclpy.node import Node
from geometry_msgs.msg import Point, PoseArray, Pose, PoseStamped
from nav_msgs.msg import OccupancyGrid, Path
from visualization_msgs.msg import Marker, MarkerArray
from std_msgs.msg import Header
from builtin_interfaces.msg import Duration


# ─────────────────────────────────────────────────────────────────────────────
# Pure-Python A* on a 2-D numpy grid
# ─────────────────────────────────────────────────────────────────────────────

def astar_2d(
    grid: np.ndarray,
    start: tuple[int, int],
    goal: tuple[int, int],
    obstacle_threshold: int = 50,
) -> Optional[list[tuple[int, int]]]:
    """
    A* on a 2-D occupancy array.

    Args:
        grid: 2-D int array, values 0-100 (100 = occupied, -1 = unknown).
        start, goal: (row, col) integer grid coordinates.
        obstacle_threshold: cells with value >= this are impassable.

    Returns:
        List of (row, col) from start to goal, or None if no path exists.
    """
    rows, cols = grid.shape

    def in_bounds(r, c):
        return 0 <= r < rows and 0 <= c < cols

    def passable(r, c):
        v = grid[r, c]
        return v != -1 and v < obstacle_threshold

    def h(a, b):
        return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2)

    # 8-connected neighbourhood
    DIRS = [(-1, 0), (1, 0), (0, -1), (0, 1),
            (-1, -1), (-1, 1), (1, -1), (1, 1)]
    COSTS = [1.0, 1.0, 1.0, 1.0,
             math.sqrt(2), math.sqrt(2), math.sqrt(2), math.sqrt(2)]

    open_heap: list = []  # (f, (r, c))
    heapq.heappush(open_heap, (h(start, goal), start))
    came_from: dict[tuple, Optional[tuple]] = {start: None}
    g: dict[tuple, float] = {start: 0.0}

    while open_heap:
        _, current = heapq.heappop(open_heap)

        if current == goal:
            path = []
            node = goal
            while node is not None:
                path.append(node)
                node = came_from[node]
            return path[::-1]

        for (dr, dc), cost in zip(DIRS, COSTS):
            nb = (current[0] + dr, current[1] + dc)
            if not in_bounds(*nb) or not passable(*nb):
                continue
            ng = g[current] + cost
            if nb not in g or ng < g[nb]:
                g[nb] = ng
                f = ng + h(nb, goal)
                came_from[nb] = current
                heapq.heappush(open_heap, (f, nb))

    return None  # no path found


# ─────────────────────────────────────────────────────────────────────────────
# Bresenham line-of-sight check
# ─────────────────────────────────────────────────────────────────────────────

def bresenham_los(
    grid: np.ndarray,
    p1: tuple[int, int],
    p2: tuple[int, int],
    obstacle_threshold: int = 50,
) -> bool:
    """Return True if the straight line p1→p2 is obstacle-free."""
    r0, c0 = p1
    r1, c1 = p2
    dr = abs(r1 - r0)
    dc = abs(c1 - c0)
    sr = 1 if r1 > r0 else -1
    sc = 1 if c1 > c0 else -1
    r, c = r0, c0
    err = dr - dc

    while True:
        v = grid[r, c]
        if v == -1 or v >= obstacle_threshold:
            return False
        if r == r1 and c == c1:
            return True
        e2 = 2 * err
        if e2 > -dc:
            err -= dc
            r += sr
        if e2 < dr:
            err += dr
            c += sc


# ─────────────────────────────────────────────────────────────────────────────
# Visibility string-pulling
# ─────────────────────────────────────────────────────────────────────────────

def string_pull(
    grid: np.ndarray,
    path: list[tuple[int, int]],
) -> list[tuple[int, int]]:
    """
    Reduce an A* path to the minimum set of waypoints that preserves
    obstacle-free LoS between consecutive nodes (Theta*/string-pulling).

    This gives us the "corner anchors" where relay drones should station
    themselves: they are exactly the points where the path bends around
    a building corner.
    """
    if len(path) < 3:
        return path

    anchors = [path[0]]
    i = 0
    while i < len(path) - 1:
        # Greedy backward search: try to connect current anchor to the
        # furthest visible point in the remaining path.
        j = len(path) - 1
        while j > i + 1:
            if bresenham_los(grid, path[i], path[j]):
                break
            j -= 1
        anchors.append(path[j])
        i = j

    return anchors


# ─────────────────────────────────────────────────────────────────────────────
# Comm-range enforcement
# ─────────────────────────────────────────────────────────────────────────────

def enforce_comm_range(
    waypoints_world: list[np.ndarray],
    comm_range: float,
) -> list[np.ndarray]:
    """
    If two consecutive waypoints are farther than comm_range apart, insert
    equally-spaced intermediate points along the segment so every adjacent
    pair is within comm_range.
    """
    if len(waypoints_world) < 2:
        return waypoints_world

    result = [waypoints_world[0]]
    for i in range(1, len(waypoints_world)):
        p1 = result[-1]
        p2 = waypoints_world[i]
        dist = float(np.linalg.norm(p2 - p1))
        if dist > comm_range:
            n_segs = math.ceil(dist / comm_range)
            for k in range(1, n_segs):
                t = k / n_segs
                result.append(p1 + t * (p2 - p1))
        result.append(p2)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Grid inflation
# ─────────────────────────────────────────────────────────────────────────────

def inflate_obstacles(
    grid: np.ndarray,
    radius_cells: int,
    obstacle_threshold: int = 50,
) -> np.ndarray:
    """
    Expand every occupied cell by radius_cells in all directions.
    Uses a simple square-kernel dilation without scipy dependency.
    """
    occupied = (grid >= obstacle_threshold) | (grid == -1)
    inflated = np.zeros_like(occupied)
    rows, cols = grid.shape
    for dr in range(-radius_cells, radius_cells + 1):
        for dc in range(-radius_cells, radius_cells + 1):
            shifted = np.zeros_like(occupied)
            r0, r1 = max(0, -dr), min(rows, rows - dr)
            c0, c1 = max(0, -dc), min(cols, cols - dc)
            sr0, sr1 = max(0, dr), min(rows, rows + dr)
            sc0, sc1 = max(0, dc), min(cols, cols + dc)
            shifted[sr0:sr1, sc0:sc1] = occupied[r0:r1, c0:c1]
            inflated |= shifted

    result = grid.copy()
    result[inflated] = 100
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Node
# ─────────────────────────────────────────────────────────────────────────────

class PathPlanningNode(Node):

    FALLBACK_SPACING = 15.0  # metres between relay points when no map

    def __init__(self):
        super().__init__('path_planning_node')

        # ── Parameters ───────────────────────────────────────────────
        self.declare_parameter('flight_height', 1.5)
        self.declare_parameter('safety_margin', 0.5)   # metres
        self.declare_parameter('comm_range',    70.0)  # metres
        self.declare_parameter('replan_dist',   2.0)   # metres
        self.declare_parameter('base_x',        0.0)
        self.declare_parameter('base_y',        0.0)

        self.flight_height = self.get_parameter('flight_height').value
        self.safety_margin = self.get_parameter('safety_margin').value
        self.comm_range    = self.get_parameter('comm_range').value
        self.replan_dist   = self.get_parameter('replan_dist').value
        self.base          = np.array([
            self.get_parameter('base_x').value,
            self.get_parameter('base_y').value,
        ])

        # ── State ────────────────────────────────────────────────────
        self.grid_map: Optional[OccupancyGrid] = None
        self.agv_pos:  Optional[np.ndarray]    = None
        self.last_agv_replan: Optional[np.ndarray] = None
        self.relay_waypoints: list[np.ndarray] = []

        # ── Subscriptions ────────────────────────────────────────────
        self.create_subscription(
            OccupancyGrid, '/projected_map',
            self._map_cb, 1)
        self.create_subscription(
            Point, '/AGV/pose',
            self._agv_cb, 10)

        # ── Publishers ───────────────────────────────────────────────
        self.wp_pub   = self.create_publisher(PoseArray,   '/path_planning/relay_waypoints', 10)
        self.path_pub = self.create_publisher(Path,        '/path_planning/agv_path',        10)
        self.viz_pub  = self.create_publisher(MarkerArray, '/path_planning/viz',              10)

        # ── Timer ────────────────────────────────────────────────────
        self.create_timer(0.5, self._tick)

        self.get_logger().info(
            f"PathPlanningNode ready  |  safety_margin={self.safety_margin} m  "
            f"|  comm_range={self.comm_range} m"
        )

    # ── Callbacks ─────────────────────────────────────────────────────────

    def _map_cb(self, msg: OccupancyGrid):
        self.grid_map = msg
        self.last_agv_replan = None  # force replan with new map

    def _agv_cb(self, msg: Point):
        self.agv_pos = np.array([msg.x, msg.y])

    # ── Main tick ─────────────────────────────────────────────────────────

    def _tick(self):
        if self.agv_pos is None:
            return

        # Replan only when AGV has moved enough
        if (self.last_agv_replan is not None
                and np.linalg.norm(self.agv_pos - self.last_agv_replan)
                < self.replan_dist):
            return

        self.last_agv_replan = self.agv_pos.copy()

        if self.grid_map is not None:
            waypoints = self._plan_with_map()
        else:
            waypoints = self._plan_fallback()

        self.relay_waypoints = waypoints
        self._publish(waypoints)

    # ── Map-based planning ────────────────────────────────────────────────

    def _plan_with_map(self) -> list[np.ndarray]:
        """
        Full A* planning on the 2-D projected occupancy grid.

        Steps:
          1. Convert OccupancyGrid to numpy array.
          2. Inflate obstacles by safety_margin.
          3. Convert base and AGV positions to grid coordinates.
          4. Run A*.
          5. String-pull to extract corner anchor waypoints.
          6. Enforce comm_range between consecutive anchors.
          7. Convert back to world coordinates.
        """
        msg = self.grid_map
        info = msg.info
        res  = info.resolution
        w, h = info.width, info.height
        ox   = info.origin.position.x
        oy   = info.origin.position.y

        # 1. OccupancyGrid → 2-D array (row = y, col = x)
        raw = np.array(msg.data, dtype=np.int16).reshape((h, w))

        # 2. Inflate obstacles
        rad_cells = max(1, int(math.ceil(self.safety_margin / res)))
        grid = inflate_obstacles(raw, rad_cells)

        # 3. World → grid coords
        def world_to_grid(wx, wy):
            col = int((wx - ox) / res)
            row = int((wy - oy) / res)
            col = max(0, min(w - 1, col))
            row = max(0, min(h - 1, row))
            return row, col

        def grid_to_world(row, col):
            wx = ox + (col + 0.5) * res
            wy = oy + (row + 0.5) * res
            return np.array([wx, wy])

        start_rc = world_to_grid(self.base[0],    self.base[1])
        goal_rc  = world_to_grid(self.agv_pos[0], self.agv_pos[1])

        # 4. A*
        path_rc = astar_2d(grid, start_rc, goal_rc)
        if path_rc is None:
            self.get_logger().warn("A* found no path — falling back to straight line.")
            return self._plan_fallback()

        # 5. String-pull → corner anchors
        anchors_rc = string_pull(grid, path_rc)

        # 6. Convert to world; drop first (base) and last (AGV) node
        anchors_world = [grid_to_world(*rc) for rc in anchors_rc]
        relay_world   = anchors_world[1:-1]  # interior corners only

        # 7. Enforce comm_range (add midpoints if gap > comm_range)
        all_world = [self.base] + relay_world + [self.agv_pos]
        all_world = enforce_comm_range(all_world, self.comm_range)
        relay_world = all_world[1:-1]

        # Publish full path for RViz
        self._publish_path([grid_to_world(*rc) for rc in path_rc])

        return relay_world

    # ── Fallback planning (no map) ────────────────────────────────────────

    def _plan_fallback(self) -> list[np.ndarray]:
        """
        When no OccupancyGrid is available (e.g., empty world), distribute
        relay points uniformly along the straight line from base to AGV.
        """
        vec  = self.agv_pos - self.base
        dist = float(np.linalg.norm(vec))

        if dist < 1e-3:
            return []

        n_relays = max(0, int(dist / self.comm_range))
        waypoints = []
        for i in range(1, n_relays + 1):
            t = i / (n_relays + 1)
            waypoints.append(self.base + t * vec)

        return waypoints

    # ── Publishers ────────────────────────────────────────────────────────

    def _publish(self, waypoints: list[np.ndarray]):
        """Publish relay waypoints as PoseArray and RViz markers."""
        now = self.get_clock().now().to_msg()

        # PoseArray
        pa = PoseArray()
        pa.header.frame_id = 'world'
        pa.header.stamp    = now
        for wp in waypoints:
            p = Pose()
            p.position.x = float(wp[0])
            p.position.y = float(wp[1])
            p.position.z = self.flight_height
            p.orientation.w = 1.0
            pa.poses.append(p)
        self.wp_pub.publish(pa)

        # RViz markers
        ma = MarkerArray()
        delete = Marker()
        delete.action = Marker.DELETEALL
        delete.header.frame_id = 'world'
        delete.header.stamp    = now
        ma.markers.append(delete)

        for i, wp in enumerate(waypoints):
            m = Marker()
            m.header.frame_id = 'world'
            m.header.stamp    = now
            m.ns   = 'relay_anchors'
            m.id   = i
            m.type = Marker.CYLINDER
            m.action = Marker.ADD
            m.pose.position.x = float(wp[0])
            m.pose.position.y = float(wp[1])
            m.pose.position.z = self.flight_height
            m.pose.orientation.w = 1.0
            m.scale.x = m.scale.y = 0.5
            m.scale.z = self.flight_height * 2
            m.color.r = 0.2
            m.color.g = 0.8
            m.color.b = 0.2
            m.color.a = 0.6
            m.lifetime = Duration(sec=2)
            ma.markers.append(m)

        self.viz_pub.publish(ma)

        self.get_logger().info(
            f"Published {len(waypoints)} relay anchor(s)  |  "
            f"AGV=({self.agv_pos[0]:.1f},{self.agv_pos[1]:.1f})",
            throttle_duration_sec=2.0
        )

    def _publish_path(self, path_world: list[np.ndarray]):
        """Publish the raw A* path as nav_msgs/Path for RViz."""
        msg = Path()
        msg.header.frame_id = 'world'
        msg.header.stamp    = self.get_clock().now().to_msg()
        for wp in path_world:
            ps = PoseStamped()
            ps.header = msg.header
            ps.pose.position.x = float(wp[0])
            ps.pose.position.y = float(wp[1])
            ps.pose.position.z = self.flight_height
            ps.pose.orientation.w = 1.0
            msg.poses.append(ps)
        self.path_pub.publish(msg)


# ─────────────────────────────────────────────────────────────────────────────

def main(args=None):
    rclpy.init(args=args)
    node = PathPlanningNode()
    try:
        rclpy.spin(node)
    except KeyboardInterrupt:
        pass
    finally:
        node.destroy_node()
        rclpy.shutdown()


if __name__ == '__main__':
    main()
