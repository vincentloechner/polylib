import pypolylib_core as pl
from pypolylib import LBLRead
pl.LBLRead = LBLRead


# This is a test python program (inspired from ZAlloc3d.in: IMPACT paper example)

# Create an LBL from a string:
a = pl.LBLRead("{(i+j)  | 0 <= i <= 10, i = 2j}")

# should be (k is arbitrary and projected out, but I added constraint k>=0):
# a = pl.LBLRead("{(i+j)  | 0 <= i <= 10, i = 2j, k >= 0}")

print("a =", a)

# should be something like:
# {(3i)  | 0 <= i <= 5}
