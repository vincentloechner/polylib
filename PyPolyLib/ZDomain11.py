import pypolylib_core as pl
from pypolylib import LBLRead
pl.LBLRead = LBLRead


# This is a test python program (inspired from ZDomain11.in)

# Create an LBL from a string:
a = pl.LBLRead("{(2i, 0, 52i+66j) | 198k <= -150i+j, 264k >= -200i+j, 9j <= 2i, 5j >= i, i <= 1000}")

print("a =", a)
# should be something like:
# {(2i, 0, 52i+66j) | -150i-198j+k >= 0, 200i+264j-k >= 0, 2i-9k>=0, i<=1000, i<=5k} union

# Eliminate existential variable k:

b = a.zdomain()
print("a.zdomain =", b)
# should be something like:
# {2i, 0, 52i+66j) | i <= 1000, -674i-891j>=0, 333i+440j>=0}

# Check inclusion:

# check if a is included in b
check1 = a.included(b)
# check if b is included in a
check2 = b.included(a)
# or: check1 = pl.LBLincluded(a, b) ...

print("Check inclusion (both ways):", check1, check2)
