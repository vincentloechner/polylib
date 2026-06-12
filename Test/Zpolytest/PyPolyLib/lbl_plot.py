"""
lbl_plot.py — Visualisation de Z-polyèdres avec matplotlib.

Utilise plot_polyhedron.py du prof pour afficher les LBL 2D.
"""

import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__), '../Verif'))

from plot_polyhedron import plot_polyhedra
import pypolylib_core as pl


def lbl_plot(lbl_obj):
    """
    Affiche un LBL 2D avec matplotlib.

    Extrait les contraintes du polyèdre et appelle plot_polyhedra.
    Fonctionne uniquement pour les LBL de dimension 2.

    Args:
        lbl_obj: Objet LBL Python (classe LBL de pypolylib.py)
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
            print(f"Attention : dimension {dim}, seule la dimension 2 est supportée.")
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
                # Égalité : aussi ajouter -ax - by - c >= 0
                inequalities.append((-a, -b, -c))

        list_of_polyhedra.append(inequalities)
        node = node.next

    if not list_of_polyhedra:
        print("Rien à afficher.")
        return

    plot_polyhedra(list_of_polyhedra)
