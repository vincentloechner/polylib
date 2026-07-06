# basic pypolylib tests

from pypolylib import LBL, Transfo

"""--------------------- Universe/empty LBL ----------------------"""
univ = LBL("{(x,y)|}")
assert str(univ) == "{(i, j) | }"
empty = LBL("{(x,y)|1 < 0}")
assert str(empty) == "{(i, j) | <empty>}"

### ERROR : empty LBL has not necessarily an HNF lat !!!

"""--------------------- Basic LBL manipulation ---------------------"""
a = LBL("{(2i, 2i+j) | 1 <= i <= 5, 1 <= j <= i}")
b = LBL("{(i,j) | j < 5}")
assert str(a) == "{(2i, j) | i <= 5, -2i+j - 1 >= 0, 3i-j >= 0}"
assert str(b) == "{(i, j) | j <= 4}"

# a is part of the universe, and empty is part of a:
assert a in univ
assert empty in a

# union
assert str(a + b) == """{(2i, j) | i <= 5, -2i+j - 1 >= 0, 3i-j >= 0} UNION
  {(i, j) | j <= 4}"""
assert a+b == b.union(a) == a.union(b)

# difference
assert str(a - b) == "{(2i, j) | i <= 5, -2i+j - 1 >= 0, 3i-j >= 0, j >= 5}"
assert a-b == a.difference(b)

# intersection
assert str(a * b) == "{(2i, j) | -2i+j - 1 >= 0, 3i-j >= 0, j <= 4}"
assert a*b == a.intersection(b)

# Z-domain
assert a.zdomain() == a

"""--------------------- Iterate over LBLs ---------------------"""
sA = set(a)
assert len(sA) == 15
for z in a:
  assert z in sA

"""--------------------- Errors over LBLs ---------------------"""
# incompatible dimensions
c = LBL("{(i, j, k) | 0 <= i <=5, 0 <= j <= i, 0 <= k <= i}")  
try: a+c; assert False
except ValueError: pass

# enumerate unbounded LBL
try: list(b); assert False
except ValueError: pass

# parse error
try: d = LBL("{nothing}"); assert False
except ValueError: pass

# parse error (UTF-8)
try: d = LBL("{(x)|x≤2}"); assert False
except ValueError: pass

"""--------------------- Transformations ---------------------"""
t = Transfo("(i,j -> 3i+1, i+j)")
assert str(t) == "(i, j -> 3i+1, i+j)"
c = a.image(t)      # c = image of a under t
d = c.preimage(t)   # preimage of c under f
assert d == a

print("all tests passed")
