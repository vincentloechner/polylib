"""
lbl_plot.py — Visualizing Z-polyhedra with matplotlib.
Typical usage:
    from pypolylib import LBL
    a = LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
    a.plot()
"""
import pypolylib_core as pl
# from plot_polyhedron import plot_polyhedra
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from shapely.geometry import Point, MultiPoint


def lbl_plot(lbl_obj, show_points=True):
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
        print("Nothing to display.")
        return
    # Get dimension from first node
    dim = node.Lat.nbrows - 1

    if dim == 2:
        _plot_2d(lbl_obj, show_points)
    elif dim == 3:
        _plot_3d(lbl_obj, show_points)
    else:
        print(f"Warning: dimension {dim} is not supported (only 2D and 3D).")


def _plot_2d(lbl_obj, show_points):
    """2D plot using plot_polyhedra (polygon + integer points)."""

    nb_nodes = 0
    node = lbl_obj._lbl
    while node is not None:
        nb_nodes += 1
        node = node.next

    ax = plt.subplots()[1]
    colors = plt.cm.tab10(np.linspace(0, 1, nb_nodes))

    node_num = 0
    node = lbl_obj._lbl
    while node is not None:
        p = node.P
        # loop on polyhedra inside node
        while p is not None:

            # need to compute the image of poly by Lat!
            poly = pl.PolyhedronImage(node.Lat, p)
            if poly.nbbid != 0:
                print("Warning: trying to plot an unbounded polyhedron. Canceling.")
                return

            vertices = []
            rays = poly.rays
            for v in range(rays.nbrows):
                # no need, this is always 1 (no bid lines):
                # v_type = rays[v, 0]
                a = rays[v, 1]
                b = rays[v, 2]
                c = rays[v, 3]
                if c == 0:
                    print("Warning: trying to plot an unbounded polyhedron. Canceling.")
                    return
                vertices.append((a/c, b/c))
            plot_convex(vertices, ax, colors[node_num])
            p = p.next
        node = node.next
        node_num += 1
    if show_points:
        points = list(lbl_obj)
        xs, ys = zip(*points)
        ax.scatter(xs, ys, s=20) #, color=colors[k])

    ax.set_aspect('equal')
    ax.grid(True, linestyle='--', alpha=0.4)
    plt.show()

def _plot_3d(lbl_obj, show_points):
    """3D scatter plot of integer points (interactive: rotate/zoom)."""

    pts = list(lbl_obj)
    if not pts:
        print("Nothing to display.")
        return

    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    zs = [p[2] for p in pts]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(xs, ys, zs, c='steelblue', marker='o', s=30)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title("LBL 3D")
    plt.tight_layout()
    plt.show()


def plot_convex(vertices, ax, col):
    mp = MultiPoint(vertices)
    region = mp.convex_hull
    if region is None or region.is_empty or region.geom_type not in ("Polygon", "MultiPolygon"):
        return

    px, py = region.exterior.xy

    # Polygon
    ax.plot(px, py, color=col, linewidth=2)
    ax.fill(px, py, color=col, alpha=0.15)
