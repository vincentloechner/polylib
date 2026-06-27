"""
lbl_plot.py — Visualizing Z-polyhedra with matplotlib.
Typical usage:
    from pypolylib import LBLRead
    a = LBLRead("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
    a.plot()
"""
import sys
import os
from plot_polyhedron import plot_polyhedra
import pypolylib_core as pl


def lbl_plot(lbl_obj):
    """
    Plots a 2D or 3D LBL using matplotlib.
    - 2D: displays the polyhedron region + integer points.
    - 3D: displays integer points as a 3D scatter plot (interactive).
    Nodes of other dimensions are ignored with a warning.
    Args:
        lbl_obj (LBL): Python LBL object (from pypolylib.py)
    """
    # Detect dimension from first valid node
    dim = None
    node = lbl_obj._lbl
    while node is not None:
        if node.P is not None:
            dim = node.P.dimension
            break
        node = node.next

    if dim is None:
        print("Nothing to display.")
        return

    if dim == 2:
        _plot_2d(lbl_obj)
    elif dim == 3:
        _plot_3d(lbl_obj)
    else:
        print(f"Warning: dimension {dim} is not supported (only 2D and 3D).")


def _plot_2d(lbl_obj):
    """2D plot using plot_polyhedra (polygon + integer points)."""
    list_of_polyhedra = []
    node = lbl_obj._lbl
    while node is not None:
        poly = node.P
        # Boucle interne : union de polyèdres dans le même node
        while poly is not None:
            if poly.dimension != 2:
                poly = poly.next
                continue
            cmat = poly.constraints
            nb_constraints = poly.nbconstraints
            inequalities = []
            for r in range(nb_constraints):
                eq_type = pl.MatrixGetValue(cmat, r, 0)
                a = pl.MatrixGetValue(cmat, r, 1)
                b = pl.MatrixGetValue(cmat, r, 2)
                c = pl.MatrixGetValue(cmat, r, 3)
                inequalities.append((a, b, c))
                if eq_type == 0:
                    inequalities.append((-a, -b, -c))
            list_of_polyhedra.append(inequalities)
            poly = poly.next
        node = node.next
    if not list_of_polyhedra:
        print("Nothing to display.")
        return
    plot_polyhedra(list_of_polyhedra)


def _plot_3d(lbl_obj):
    """3D scatter plot of integer points (interactive: rotate/zoom)."""
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

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
