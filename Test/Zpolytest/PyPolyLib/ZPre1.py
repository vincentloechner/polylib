import pypolylib
from pypolylib import LBLRead, Transfo

# Preimage of a Z-polyhedron
# inspired from ZPre1.in

# LBL = (2i, j) avec domaine 1 <= i,j <= 100
a = LBLRead("{(2i, j) | 1 <= i <= 100, 1 <= j <= 100}")
print("a =", a)

# fonction (i,j) -> (2i+10, 3i+j+1)
f = Transfo("(i,j -> 2i+10, 3i+j+1)")

# preimage
b = a.preimage(f)
print("preimage =", b)
