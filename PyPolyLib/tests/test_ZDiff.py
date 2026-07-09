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
    u     = diff + inter
    assert u.included(a)
    assert a.included(u)


# list of (name:str, [a:str(lbl), b:str(lbl)])
tests = [
    ("ZDiff1",
     ["{(i, j) | 1 <= i <= 10, 1 <= j <= 10}",
      "{(i, j) | 5 <= i <= 15, 5 <= j <= 15}"]),

    ("ZDiff2",
     ["{(i, j) | 0 <= i <= 5, 0 <= j <= i}",
      "{(i, j) | 5 <= i <= 7, 0 <= j <= i}"]),

    ("ZDiff3",
     ["{(i, j) | 1 <= i <= 20, 1 <= j <= 20}",
      "{(i, j) | 10 <= i <= 30, 10 <= j <= 30}"]),

    ("ZDiffInter6a",
     ["{(i, j)  | 1 <= i <= 20, 1 <= j <= 20}",
      "{(2i, j) | 1 <= i <= 50, j = 10}"]),

    ("ZDiff_triangle",
     ["{(i, j) | 0 <= i <= 8, 0 <= j <= i}",
      "{(2i, j) | 0 <= i <= 4, 0 <= j <= 2i}"]),
]



from do_tests import do_tests
do_tests(tests, test_func)

