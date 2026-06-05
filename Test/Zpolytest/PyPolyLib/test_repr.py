import pypolylib_core as pl
from pypolylib import LBLRead

a = LBLRead("{(i+j) | 0 <= i <= 10, i = 2j}")

lbl = a._lbl
lat = lbl.Lat
poly = lbl.P

print("=== LATTICE ===")
print("nbrows:", lat.nbrows, "nbcols:", lat.nbcolumns)
for i in range(lat.nbrows):
    for j in range(lat.nbcolumns):
        print(f"  [{i}][{j}] =", pl.MatrixGetValue(lat, i, j))

print("=== POLYHEDRON ===")
print("dimension:", poly.dimension)
print("nbconstraints:", poly.nbconstraints)
for i in range(poly.nbconstraints):
    for j in range(poly.dimension + 2):
        print(f"  [{i}][{j}] =", pl.MatrixGetValue(poly.constraints, i, j))
