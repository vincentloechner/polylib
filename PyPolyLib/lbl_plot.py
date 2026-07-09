"""
lbl_plot.py — Visualizing Z-polyhedra with matplotlib.
Typical usage:
    from pypolylib import LBL
    a = LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
    a.plot()
"""
import sys
import pypolylib_core as pl

# plotting libs:
import math
import numpy as np
from colorsys import hsv_to_rgb
# # 2D -> not used anymore
# import matplotlib.pyplot as plt
# from shapely.geometry import MultiPoint
# 3D
import pyvista as pv
from scipy.spatial import ConvexHull

class WindowColorGenerator:
    def __init__(self):
        self.hue = 0.0      # for color
        self._3D = False    # is there a 3D polyhedron in this window?
        self.plot = pv.Plotter()

    def color(self):
        self.hue = (self.hue + 0.618033988749895) % 1.0
        return hsv_to_rgb(self.hue, .55, .9)

# this is not very nice, but it's the best way to do it.
_lbl_plot_window = None

def lbl_plot(lbl, *args, show_points=True, subplot=False, **kwargs):
    """
    Plots an LBL using pyvista.

    Args:
    - lbl (LBL): Python LBL object (from pypolylib.py)
    - show_points (boolean): draw points that are part of the LBL
    - subplot (boolean): set this to True if you want to plot another LBL
        in the same window (it will be rendered later in that case)
    - extra arguments will be passed to the call to plotter.show()
    """
    global _lbl_plot_window

    node = lbl._lbl
    if node is None:
        return
    if _lbl_plot_window is None:
        _lbl_plot_window = WindowColorGenerator()
    plotter = _lbl_plot_window.plot

    # loop over lattices and polyhedral domains to plot them
    while node is not None:
        p = node.P
        if node.Lat.nbrows >= 4:
            _lbl_plot_window._3D = True
        lat = _lat3D(node.Lat) # always project to 3D, one way or another!
        # loop over polyhedra inside node
        while p is not None:
            # each single LBL is displayed in a different color:
            color = _lbl_plot_window.color()
            mesh = _poly2pyvista(pl.PolyhedronImage(lat, p))
            if mesh:
                plotter.add_mesh(mesh,
                                show_edges=True,
                                opacity=.25, color=color)
            if show_points:
                points = list(lbl._iter_single_lbl(lat, p, set()))
                plotter.add_mesh(
                    pv.PolyData(np.asarray(points, dtype=float)),
                    opacity=.5, color=color,
                    point_size=16, render_points_as_spheres=True
                )
            p = p.next
        node = node.next

    # if subplot is True do not render anything yet: the next call will add
    # some LBLs over this one in the same window and do the rendering
    if not subplot:
        if not _lbl_plot_window._3D:
            # this is a 2D plot in a 3D scene, disable z-axis:
            plotter.enable_2d_style()
            plotter.view_xy()
            plotter.reset_camera()

        # set_grid_ticks_1(plotter)        
        plotter.show_grid()
        # show_my_grid(plotter)
        plotter.show(**kwargs)

        # you can not keep older windows open once you quit the main
        # graphical loop, just reopen them if needed.
        if not kwargs:
            # ensure everything is closed properly if not running interactive
            pv.close_all()
        _lbl_plot_window = None # will open a new window in a next call

# original version:
    # if dim == 2:
    #     _plot_2d(lbl, show_points, subplot)
    # elif dim == 3:
    #     _plot_3d(lbl, show_points, subplot)
    # else:
    #     print(f"Plot: dimension {dim} is not supported (only 2D and 3D).",
    #           file=sys.stderr)


# def _plot_2d(lbl_obj, show_points, subplot):
#     """2D plot (polygon + integer points)."""
#     ax = plt.subplots()[1]

#     nb_nodes = 0
#     node = lbl_obj._lbl
#     while node is not None:
#         nb_nodes += 1
#         node = node.next
#     colors = plt.cm.tab10(np.linspace(0, 1, nb_nodes))

#     node = lbl_obj._lbl
#     for node_num in range(nb_nodes):
#         p = node.P
#         # loop on polyhedra inside node
#         while p is not None:
#             vertices = _get_vertices(pl.PolyhedronImage(node.Lat, p))
#             plot_convex_2D(vertices, ax, colors[node_num])
#             p = p.next
#         node = node.next

#     # need to fix color here:
#     if show_points:
#         points = list(lbl_obj)
#         xs, ys = zip(*points)
#         ax.scatter(xs, ys, s=30) #, color=colors[k])

#     ax.set_aspect('equal')
#     ax.grid(True, linestyle='--', alpha=0.4)
#     if not subplot:
#         plt.show()


# def plot_convex_2D(vertices, ax, color):
#     region = MultiPoint(vertices).convex_hull

#     if region.is_empty:
#         return

#     if region.geom_type == "Polygon":
#         x, y = region.exterior.xy
#         ax.plot(x, y, color=color, linewidth=2)
#         ax.fill(x, y, color=color, alpha=0.15)

#     elif region.geom_type == "LineString":
#         x, y = region.xy
#         ax.plot(x, y, color=color, linewidth=2)

# test:
# def set_grid_ticks_1(plotter):
#     xmin, xmax, ymin, ymax, zmin, zmax = plotter.bounds
#     xmin = math.floor(xmin)
#     xmax = math.ceil(xmax)
#     ymin = math.floor(ymin)
#     ymax = math.ceil(ymax)
#     zmin = math.floor(zmin)
#     zmax = math.ceil(zmax)
#     cube_axes_actor = pv.CubeAxesActor(plotter.camera)
#     cube_axes_actor.n_xlabels = xmax - xmin
#     cube_axes_actor.n_ylabels = ymax - ymin
#     # cube_axes_actor.n_zlabels = zmax - zmin
#     actor, property = plotter.add_actor(cube_axes_actor)

# test:
# def show_my_grid(plotter):
#     xmin, xmax, ymin, ymax, zmin, zmax = plotter.bounds
#     xmin = math.floor(xmin)
#     xmax = math.ceil(xmax)
#     ymin = math.floor(ymin)
#     ymax = math.ceil(ymax)
#     zmin = math.floor(zmin)
#     zmax = math.ceil(zmax)
#     for x in range(xmin, xmax + 1):
#         plotter.add_lines(np.array([[x, ymin, 0], [x, ymax, 0]]), color="lightgray")

#     for y in range(ymin, ymax + 1):
#         plotter.add_lines(np.array([[xmin, y, 0], [xmax, y, 0]]), color="lightgray")

#     for z in range(zmin, zmax + 1):
#         plotter.add_lines(np.array([[xmin, ymin, z], [xmax, ymin, z]]), color="lightgray")


def _poly2pyvista(poly):
    """Transform a PolyLib 3D Polyhedron into a PyVista polyhedron."""

    # vertices: as the list of tuples of FP coordinates
    vertices = _get_vertices(poly)
    vertices = np.asarray(vertices, dtype=float)

    # faces: as the (flat) list of faces [num_vertices, vertex0, vertex1, ...]
    faces = []
    dim = poly.dimension
    constraint = poly.constraint
    ray = poly.ray
    # scan constraints to build faces:
    for c in range(poly.nbconstraints):
        normal = tuple(constraint[c, d+1] for d in range(dim))
        face = []
        for r in range(poly.nbrays):
            sat = sum(constraint[c, d+1]*ray[r, d+1] for d in range(dim+1))
            # does ray r saturate constraint c?
            if sat == 0:
                face.append(r)
        if len(face) >= 3:
            face = _order_face_indices(vertices, face, normal)
            faces.extend([len(face)] + face)
    if faces:
        return pv.PolyData(vertices, faces)


def _order_face_indices(vertices, face_indices, normal):
    face_indices = np.asarray(face_indices, dtype=int)
    ordered_local = _order_face(vertices[face_indices], normal)
    return [face_indices[i] for i in ordered_local]


def _order_face(vertices, normal):
    """Order coplanar 3D vertices cyclically."""
    normal = np.asarray(normal, dtype=float)
    # Normalize normal
    normal /= np.linalg.norm(normal)
    if len(vertices) == 3:
        order = np.array([0,1,2])
    else:
        # Build an orthonormal basis (u,v,n): choose a vector not parallel to n
        if abs(normal[0]) < 0.9:
            tmp = np.array([1., 0., 0.])
        else:
            tmp = np.array([0., 1., 0.])
        u = np.cross(normal, tmp)
        u /= np.linalg.norm(u)
        v = np.cross(normal, u)

        # Project points onto the face plane
        center = vertices.mean(axis=0)
        centered = vertices - center
        pts2d = np.column_stack((centered @ u, centered @ v))

        # Convex hull in 2D
        hull = ConvexHull(pts2d)
        order = hull.vertices

    # Fix orientation
    p0, p1, p2 = vertices[order[:3]]
    face_normal = np.cross(p1 - p0, p2 - p0)
    if np.dot(face_normal, normal) < 0:
        order = order[::-1]
    return order

def _get_vertices(poly):
    """ Return the vertices of poly as a list of tuples of floats.

    return [] if poly is not bounded (and print a message to stderr).
    """

    if poly.nbbid != 0:
        print("Warning: trying to plot an unbounded polyhedron", file=sys.stderr)
        return []

    vertices = []
    dim = poly.dimension
    ray = poly.ray
    for v in range(poly.nbrays):
        # v_type = ray[v, 0] # alway 1 (no bid. ray)
        v_div = ray[v, dim + 1]
        if v_div == 0:
            print("Warning: trying to plot an unbounded polyhedron", file=sys.stderr)
            return []
        vertices.append(tuple(ray[v, x+1]/v_div for x in range(dim)))

    return vertices

def _lat3D(node_lat):
    """Change a Lattice to be a 3D representation/projection.
    
    - if lower than 3D -> set the plane z=0 (and y=0 if 1D)
    - if greater than 3D -> project out extra dimensions (TODO: could be improved)
    """
    if node_lat.nbrows == 4:
        # 3D is ok
        return node_lat
    else:
        # extend node_lat to be a 3D projection!
        lat = pl.MatrixAlloc(4, node_lat.nbcolumns)
        for i in range(min(3, node_lat.nbrows - 1)):
            for j in range(node_lat.nbcolumns):
                lat[i, j] = node_lat[i, j]
        lat[3, lat.nbcolumns - 1] = 1
        return lat
