# basic pypolylib tests

from pypolylib import LBL, Transfo, PolylibClose

def main():
  """--------------------- Universe/empty LBL ----------------------"""
  print(f"Universe/empty LBL     ", end="", flush=True)

  univ = LBL("{(x,y)|}")
  assert str(univ) == "{(i, j) | }"
  empty = LBL("{(x,y)|1 < 0}")
  assert str(empty) == "{(_, _) | <empty>}"
  print("OK")

  """--------------------- Basic LBL manipulation ---------------------"""
  print(f"2D LBL i/o             ", end="", flush=True)
  a = LBL("{(2i, 2i+j) | 1 <= i <= 5, 1 <= j <= i}")
  b = LBL("{(i,j) | j < 5}")
  assert str(a) == "{(2i, j) | i <= 5, -2i+j - 1 >= 0, 3i-j >= 0}"
  assert str(b) == "{(i, j) | j <= 4}"

  # a is part of the universe, and empty is part of a:
  assert a in univ
  assert empty in a

  print("OK")

  # operations
  print(f"Operations on 2D LBL   ", end="", flush=True)
  # union
  u = a + b
  assert str(u) == """{(2i, j) | i <= 5, -2i+j - 1 >= 0, 3i-j >= 0} UNION
  {(i, j) | j <= 4}"""
  assert u == b.union(a) == a.union(b)

  # difference
  d = a - b
  assert str(d) == "{(2i, j) | i <= 5, -2i+j - 1 >= 0, 3i-j >= 0, j >= 5}"
  assert d == a.difference(b)

  # intersection
  i = a * b
  assert str(i) == "{(2i, j) | -2i+j - 1 >= 0, 3i-j >= 0, j <= 4}"
  assert i == a.intersection(b)

  # diff * inter = {0}; diff + inter = a
  assert d * i == empty
  assert d + i == a

  # Z-domain
  assert a.zdomain() == a

  # single contains point test
  assert (2, 3) in a

  print("OK")
  """--------------------- Iterate over LBLs ---------------------"""
  print(f"Enumerate LBL          ", end="", flush=True)
  sA = set(a)             # build the (python-)set of points in a
  assert len(sA) == 15
  for z in a:
    assert z in sA        # python set inclusion test
    assert z in a         # LBL contains point test

  print("OK")
  """--------------------- Transformations ---------------------"""
  print(f"LBL transformations    ", end="", flush=True)
  t = Transfo("(i, j -> 3i+1, i+j)")       # t is a bijection
  assert str(t) == "(i, j -> 3i+1, i+j)"
  im_a = a.image(t)                       # im_a = image of a under t
  pre = im_a.preimage(t)                  # preimage of im_a under f
  assert pre == a

  t2 = Transfo("(i, j -> j, i)")
  assert str(t2 * t) == "(i, j -> i+j, 3i+1)"

  t_inv = t.inverse()
  assert str(t_inv) == "(i, j -> i-1, -i+3j+1)"
  assert t(a) == t_inv.preimage(a)

  print("OK")
  """--------------------- Value Errors over LBLs ---------------------"""
  print(f"Value errors over LBLs ", end="", flush=True)
  # incompatible dimensions
  c = LBL("{(i, j, k) | 0 <= i <=5, 0 <= j <= i, 0 <= k <= i}")
  try: a+c; assert False
  except ValueError: pass

  try: (1,2) in c; assert False
  except ValueError: pass

  # enumerate unbounded LBL
  try: list(b); assert False
  except ValueError: pass

  # parse error
  try: d = LBL("{g4r|o4ge}"); assert False
  except ValueError: pass

  # parse error (UTF-8)
  try: d = LBL("{(x)|x≤2}"); assert False
  except ValueError: pass

  print("OK")
  """--------------------- Type Errors over LBLs ---------------------"""
  print(f"Type errors over LBLs  ", end="", flush=True)
  # image of an LBL by an LBL
  try:   a.image(b); assert False
  except TypeError: pass

  # union of an LBL and a Transfo
  try:   a+t; assert False
  except TypeError: pass

  # inclusion test, LBL in Transfo:
  try:   a.included(t); assert False
  except TypeError: pass

  # contains point test, not an iterable/string:
  try:   1 in a; assert False
  except TypeError: pass
  try:   "1" in a; assert False
  except TypeError: pass

  print("OK")
  """--------------------- All done ---------------------"""
  print("all tests passed")


################################## MAIN ######################################
main()
# cleanup memory:
PolylibClose()
