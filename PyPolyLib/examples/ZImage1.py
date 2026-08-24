import pypolylib
from pypolylib import LBLRead, Transfo

# Image of a Z-polyhedron by a full-row affine function
# inspired from ZImage1.in

# LBL = lattice (i,j->2i, 2i+j) with domain {1 <= i,j <= 100}
a = LBLRead("{(2i, 2i+j) | 1 <= i <= 100, 1 <= j <= 100}")
print("a =", a)

# affine function (i,j -> 3i+1, 2i+5j)
f = Transfo("(i,j -> 3i+1, 2i+5j)")

# compute image
b = a.image(f)
print("image =", b)