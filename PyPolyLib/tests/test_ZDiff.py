"""
test_ZDiff.py — Verify ((a-b) + (a∩b) == a) for several couples of LBLs (a, b)

Usage:
    python3 test_ZDiff.py
"""
from pypolylib import LBL


# Check if (a-b) + (a∩b) == a
def test_func(a, b):
  inter = a * b
  diff  = a - b
  assert list(diff * inter) == [] # empty intersection
  u     = diff + inter
  assert u.included(a)
  assert a.included(u)


# list of (name:str, a:lbl, b:lbl)
tests = [
  ("ZDiff1",
    LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}"),
    LBL("{(i, j) | 5 <= i <= 15, 5 <= j <= 15}")),

  ("ZDiff2",
    LBL("{(i, j) | 0 <= i <= 5, 0 <= j <= i}"),
    LBL("{(i, j) | 5 <= i <= 7, 0 <= j <= i}")),

  ("ZDiff3",
    LBL("{(i, j) | 1 <= i <= 20, 1 <= j <= 20}"),
    LBL("{(i, j) | 10 <= i <= 30, 10 <= j <= 30}")),

  ("ZDiffInter6a",
    LBL("{(i, j)  | 1 <= i <= 20, 1 <= j <= 20}"),
    LBL("{(2i, j) | 1 <= i <= 50, j = 10}")),

  ("ZDiff_triangle",
    LBL("{(i, j) | 0 <= i <= 8, 0 <= j <= i}"),
    LBL("{(2i, j) | 0 <= i <= 4, 0 <= j <= 2i}")),

  # this one has a single intersection point between the two LBLs (0, 0, 0),
  # but the first one has 1030301 points
  ("ZDiff11b",
    LBL("{(2i, 5j, 52i+27j+66k) | -50 <= i <= 50, -50 <= j <= 50, -50 <= k <= 50}"),
    LBL("{(i, j, k) | i>=0, i <= 10j, 9j <= i, j+4i <= 4k, 3k <= j+3i}")),

  # inspired from the previous one, but with a larger intersection
  # a has 1030301 points
  # a*b has 49007 points
  # a-b is a union of 6 single LBLs having 981294 points = 1030301 - 49007
  ("ZDiff11c",
    LBL("{(2i, 5j, 52i+27j+66k) | -50 <= i <= 50, -50 <= j <= 50, -50 <= k <= 50}"),
    LBL("{(i, j, k) | -200 <= i <= 200, -200 <= j <= 200, -200 <= k <= 200}")),

  # odd i's in b never intersect even i's in a
  ("ZDiff11d",
    LBL("{(2i, 5j, 52i+27j+66k) | -50 <= i <= 50, -50 <= j <= 50, -50 <= k <= 50}"),
    LBL("{(2i+1, j, k) | -200 <= i <= 200, -200 <= j <= 200, -200 <= k <= 200}")),

  # only some j's intersect (5j == 2j'+3)
  ("ZDiff11e",
    LBL("{(2i, 5j, 52i+27j+66k) | -50 <= i <= 50, -50 <= j <= 50, -50 <= k <= 50}"),
    LBL("{(2i, 2j+3, k) | -200 <= i <= 200, -200 <= j <= 200, -200 <= k <= 200}")),

  # TODO:
  # in the result of ZDiff11e, the first LBL:
  # {(2i, -250, 52i+66j-1350) | j <= 50, i >= -50, i <= 50, j >= -50}
  # is included in the last one from the union:
  # {(2i, 10j, 52i+54j+66k) | i >= -50, i <= 50, j >= -25, j <= 25, k >= -50, k <= 50}
  # this could be simplified! (or not built at all, possible?)

  # simpler example in 2D (inclusion test of the first on the last of the difference):
  ("ZDiff11e_2D",
    LBL("{(2i, 5j) | -50 <= i <= 50, -50 <= j <= 50}"),
    LBL("{(2i, 2j+3) | -200 <= i <= 200, -200 <= j <= 200}")),

  # a-b
  # -> {(2i, -250) | i <= 50, i >= -50} UNION
  #   {(2i, 250) | i <= 50, i >= -50} UNION
  #   {(2i, 10j) | i >= -50, i <= 50, j >= -25, j <= 25}
  # first = LBL("{(2i, -250) | i <= 50, i >= -50}")
  # last = LBL("{(2i, 10j) | i >= -50, i <= 50, j >= -25, j <= 25}")
  # first in last
  # -> True

  # this one is from LatDiff7, to check if the lattice difference is ok:
  ("ZDiff_LatDiff7",
    LBL("{(i, 7j, k) | -50 <= i <= 50, -50 <= j <= 50, -50 <= k <= 50}"),
    LBL("{(3i+1, i+7j+3, i+2j+4k+1) | -50 <= i <= 50, -50 <= j <= 50, -50 <= k <= 50}")),

  # this one is from LatDiff10:
  ("ZDiff_LatDiff10",
    LBL("{(i, 5, k) | 0 <= i <= 20,  0 <= k <= 100}"),
    LBL("{(2i, 5i+5, 13i+33k+27) | 0 <= i <= 10, -1 <= k <= 5}")),

]


from do_tests import do_tests
do_tests(tests, test_func)

