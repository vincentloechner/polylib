"""
test_ZDisj.py — Build the disjoint union (a-b) and check that each LBL in this
list is not included in any other.

Usage:
    python3 test_ZDisj.py
"""
from pypolylib import LBL


# Check if (a-b) + (a∩b) == a
def test_func(a, b):
  diff  = a - b
  u     = diff.disjoint()
  ul    = u.sLBL_list()
  ll = len(ul)
  for i in range(ll):
    for j in range(i):
      if (ul[i] * ul[j]):
        raise Exception(f"lbl[{i}] {ul[i]} intersects lbl[{j}] {ul[j]}")


# list of (name:str, a:lbl, b:lbl)
tests = [
  ("demo",
    LBL("{(i,j) | 0 <= i <= 6, 0 <= j <= 8}"),
    LBL("{(2i, 2i+j) | 1 <= i <= 5, 1 <= j <= i}"),
  ),

]


from do_tests import do_tests
do_tests(tests, test_func)

