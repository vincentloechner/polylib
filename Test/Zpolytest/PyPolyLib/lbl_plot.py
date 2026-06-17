"""
lbl_plot.py — Visualizing Z-polyhedra with matplotlib.

Typical usage:
    from pypolylib import LBLRead
    a = LBLRead("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
    a.plot()   # affiche le LBL avec matplotlib

    # Ou directement :
    from lbl_plot import lbl_plot
    lbl_plot(a)
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../Verif'))

from plot_polyhedron import plot_polyhedra
import pypolylib_core as pl


def lbl_plot(lbl_obj):
    """
    Plots a 2D LBL using matplotlib.

    Extracts the polyhedron constraints from each node of the LBL
    and calls `plot_polyhedra` to display them.

    Works only for 2-dimensional LBLs.
    Nodes of other dimensions are ignored, and a warning message is displayed.
    Unions are displayed in different colors.

    Args:
        lbl_obj (LBL): Python LBL object (LBL class from pypolylib.py)

    Example:
        >>> from pypolylib import LBLRead
        >>> a = LBLRead(“{(i, j) | 1 <= i <= 10, 1 <= j <= 10}”)
        >>> a.plot()

        >>> # Display a union
        >>> b = LBLRead(“{(i, j) | 5 <= i <= 15, 5 <= j <= 15}”)
        >>> (a + b).plot()

    Note:
        Requires matplotlib and shapely to be installed:
            pip install matplotlib shapely
    """
    list_of_polyhedra = []
    node = lbl_obj._lbl

    while node is not None:
        poly = node.P

        if poly is None:
            node = node.next
            continue

        dim = poly.dimension
        if dim != 2:
            print(f"Warning : dimension {dim}; only dimension 2 is supported.")
            node = node.next
            continue

        cmat = poly.constraints
        nb_constraints = poly.nbconstraints

        inequalities = []
        for r in range(nb_constraints):
            eq_type = pl.MatrixGetValue(cmat, r, 0)
            a = pl.MatrixGetValue(cmat, r, 1)
            b = pl.MatrixGetValue(cmat, r, 2)
            c = pl.MatrixGetValue(cmat, r, 3)
            # PolyLib : ax + by + c >= 0
            inequalities.append((a, b, c))
            if eq_type == 0:
                # Equality: also add -ax - by - c >= 0
                inequalities.append((-a, -b, -c))

        list_of_polyhedra.append(inequalities)
        node = node.next

    if not list_of_polyhedra:
        print("Nothing to display.")
        return

    plot_polyhedra(list_of_polyhedra)