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
import numpy as np
# 2D
import matplotlib.pyplot as plt
from shapely.geometry import MultiPoint
# 3D
import pyvista as pv


def lbl_plot(lbl_obj, show_points=True, subplot=False):
    """
    Plots a 2D or 3D LBL using matplotlib.
    - 2D: displays the polyhedron region + points.
    - 3D: displays integer points as a 3D scatter plot (interactive).
    Nodes of other dimensions are ignored with a warning.
    Args:
        lbl_obj (LBL): Python LBL object (from pypolylib.py)
        show_points (boolean): draw points that are part of the LBL
    """
    node = lbl_obj._lbl
    if node is None:
        return

    _plot_3d(lbl_obj, show_points, subplot)
    # if dim == 2:
    #     _plot_2d(lbl_obj, show_points, subplot)
    # elif dim == 3:
    #     _plot_3d(lbl_obj, show_points, subplot)
    # else:
    #     print(f"Plot: dimension {dim} is not supported (only 2D and 3D).",
    #           file=sys.stderr)


def _plot_2d(lbl_obj, show_points, subplot):
    """2D plot (polygon + integer points)."""
    ax = plt.subplots()[1]

    nb_nodes = 0
    node = lbl_obj._lbl
    while node is not None:
        nb_nodes += 1
        node = node.next
    colors = plt.cm.tab10(np.linspace(0, 1, nb_nodes))

    node = lbl_obj._lbl
    for node_num in range(nb_nodes):
        p = node.P
        # loop on polyhedra inside node
        while p is not None:
            vertices = get_vertices(pl.PolyhedronImage(node.Lat, p))
            plot_convex_2D(vertices, ax, colors[node_num])
            p = p.next
        node = node.next

    # need to fix color here:
    if show_points:
        points = list(lbl_obj)
        xs, ys = zip(*points)
        ax.scatter(xs, ys, s=30) #, color=colors[k])

    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.4)
    if not subplot:
        plt.show()


def plot_convex_2D(vertices, ax, color):
    region = MultiPoint(vertices).convex_hull

    if region.is_empty:
        return

    if region.geom_type == "Polygon":
        x, y = region.exterior.xy
        ax.plot(x, y, color=color, linewidth=2)
        ax.fill(x, y, color=color, alpha=0.15)

    elif region.geom_type == "LineString":
        x, y = region.xy
        ax.plot(x, y, color=color, linewidth=2)

# this is not very nice, but best way to do it...
my_plot = None

# TODO: FIW COLORS
def _plot_3d(lbl_obj, show_points, subplot):
    """3D plot of an LBL."""
    global my_plot
    if my_plot is None:
        my_plot = pv.Plotter()
    plotter = my_plot

    nb_nodes = 0
    node = lbl_obj._lbl
    while node is not None:
        nb_nodes += 1
        node = node.next
    colors = plt.cm.tab10(np.linspace(0, 1, nb_nodes))

    node = lbl_obj._lbl
    for node_num in range(nb_nodes):
        p = node.P

        if node.Lat.nbrows == 4:
            lat = node.Lat
        else:
            # extend node.Lat to be 3D!
            lat = pl.MatrixAlloc(4, node.Lat.nbcolumns)
            for i in range(min(3, node.Lat.nbrows - 1)):
                for j in range(node.Lat.nbcolumns):
                    lat[i, j] = node.Lat[i, j]
            lat[3, lat.nbcolumns - 1] = 1

        # loop on polyhedra inside node
        while p is not None:
            mesh = poly2pyvista(pl.PolyhedronImage(lat, p))
            plotter.add_mesh(mesh,
                             show_edges=True,
                             opacity=.25, color=colors[node_num])
            if show_points:
                points = list(lbl_obj._iter_single_lbl(lat, p, set()))
                plotter.add_mesh(
                    pv.PolyData(np.asarray(points, dtype=float)),
                    opacity=.5, color=colors[node_num],
                    point_size=16, render_points_as_spheres=True
                )

            p = p.next
        node = node.next

    if not subplot:
        if lbl_obj._lbl.Lat.nbrows <= 3:
            # 2D plot in a 3D scene
            plotter.enable_2d_style()
            plotter.view_xy()
            plotter.reset_camera()
        plotter.show_grid()
        plotter.show()
        my_plot = None # will open new plot after this one is shown


def poly2pyvista(poly):

    # vertices is the list of tuples of FP coordinates
    vertices = get_vertices(poly)

    faces = []
    dim = poly.dimension
    constraint = poly.constraint
    ray = poly.ray
    for c in range(poly.nbconstraints):
        face = []
        for r in range(poly.nbrays):
            # does ray r saturate constraint c?
            sat = sum(constraint[c, d+1]*ray[r, d+1] for d in range(dim+1))
            if sat == 0:
                v_div = ray[r, dim+1]
                face.append(r)

        if len(face) >= 3:
            faces.extend([len(face)] + face)
    # faces is the list of faces as [num_vertices, vertex0, vertex1, ...]

    return pv.PolyData(vertices, faces)
    

def get_vertices(poly):    
    """ Return the vertices of poly as a list of tuples of floats.

    return [] if poly is not bounded
    """

    if poly.nbbid != 0:
        print("Warning: trying to plot an unbounded polyhedron", file=sys.stderr)
        return []

    vertices = []
    dim = poly.dimension
    ray = poly.ray
    for v in range(poly.nbrays):
        # this is always 1 (no bid lines):
        # v_type = ray[v, 0]
        v_div = ray[v, dim + 1]
        if v_div == 0:
            print("Warning: trying to plot an unbounded polyhedron", file=sys.stderr)
            return []

        vertices.append(tuple(ray[v, x+1]/v_div for x in range(dim)))
        # a = ray[v, 1]
        # b = ray[v, 2]
        # c = ray[v, 3]
        # vertices.append((a/c, b/c))

    return vertices
