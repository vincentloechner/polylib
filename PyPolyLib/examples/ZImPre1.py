import pypolylib
from pypolylib import LBLRead, Transfo

# Image then Preimage
# inspired from ZImPre1.in

# LBL = (2i, 6j) avec domaine 1 <= i,j <= 100
a = LBLRead("{(2i, 6j) | 1 <= i <= 100, 1 <= j <= 100}")
print("a =", a)

# fonction (i,j) -> (2i, 3j)
f = Transfo("(i,j -> 2i, 3j)")

# image
b = a.image(f)
print("image =", b)

# preimage
c = b.preimage(f)
print("preimage of image =", c)

# verifier inclusion
check1 = c.included(a)
check2 = a.included(c)
print("Check inclusion (both ways):", check1, check2)
