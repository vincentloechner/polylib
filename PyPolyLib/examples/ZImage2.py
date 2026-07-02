import pypolylib
from pypolylib import LBLRead, Transfo

# Image of a pointy Z-polyhedron by a projection
# inspired from ZImage2.in

# LBL = identite sur (i,j) avec domaine pointu
a = LBLRead("{(i, j) | i >= 0, i <= 10, 3j >= i, 2j <= i}")
print("a =", a)

# fonction de projection (i,j) -> i
f = Transfo("(i,j -> i)")
b = a.image(f)
print("image =", b)