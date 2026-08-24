import pypolylib
from pypolylib import LBLRead, Transfo

# Image of a Z-polyhedron by a projection from 3D to 1D
# inspired from ZImage4.in

# LBL = identite 3D avec domaine pointu
a = LBLRead("{(i, j, k) | i >= 0, i <= 100, 3j >= i, 2j <= i, 4k >= j+i, 3k <= j+i}")
print("a =", a)

# fonction (i,j,k) -> i
f = Transfo("(i,j,k -> i)")
b = a.image(f)
print("image =", b)