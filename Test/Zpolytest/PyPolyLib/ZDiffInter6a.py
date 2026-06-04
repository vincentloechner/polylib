import pypolylib
from lbl_read import LBLRead
pypolylib.LBLRead = LBLRead

# This is a test python program (inspired from ZDiffInter6a.in)

# Create two LBLs from strings:
a = pypolylib.LBLRead("{(i, j)  | 1 <= i <= 20, 1 <= j <= 20}")
b = pypolylib.LBLRead("{(2i, j) | 1 <= i <= 50, j = 10}")
# or: b = pypolylib.LBLRead("{(2i, 10) | 1 <= i <= 50}")

# in a: all integer (i,j) in [1, 20]
# in b: only even i's up to 100, and j = 10.
print("a =", a)
print("b =", b)

# inter = pypolylib.LBLintersection(a, b)
# or:
inter = a.intersection(b)
print("inter = a inter b =", inter)
# inter should be equal to: {(2i, 10) | 1 <= i <= 10}

# diff = pypolylib.LBLdifference(a, b)
# or : diff = a - b
# or:
diff = a.difference(b)
print("diff = a - b =", diff)

# diff should be equal to:
# {(i, j) | 1 <= i <= 20, 11 <= j <= 20} union
# {(i, j) | 1 <= i <= 20, 1 <= j <= 9} union
# {(1, j) | 1 <= j <= 20} union
# {(2i+1, 10) | 0 <= i <= 9}


# compute diff + inter (should be equal to original LBL, as a union)
# or : u = pypolylib.LBLunion(diff, inter)
# or : u = a + b
# or:
u = diff.union(inter)
print("u = diff union inter =", u)

# check if u is included in a
check1 = u.included(a)
# check if a is included in u
check2 = a.included(u)
# or: check1 = pypolylib.LBLincluded(u, a) ...

print("Check inclusion (both ways):", check1, check2)

