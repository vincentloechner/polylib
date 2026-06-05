import pypolylib
from pypolylib import LBLRead, Transfo

# Image of a Z-polyhedron by a projection affine function
# inspired from ZImage3.in

# LBL = lattice (2i) avec domaine 1 <= i <= 50
a = LBLRead("{(2i) | 1 <= i <= 50}")
print("a =", a)

# fonction (i) -> (i, i)
f = Transfo("(i -> i, i)")
b = a.image(f)
print("image =", b)