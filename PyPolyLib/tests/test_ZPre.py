"""
test_ZPre.py — test preimage of LBLs: image(preimage(a)) in a
"""
from pypolylib import LBL, Transfo


# test function: check that f(f^{-1}(a)) \in a
# preimage(image(a)) ⊆ a
def preimage_test_func(a, f):
    pre = a.preimage(f)
    img = pre.image(f)
    assert img in a, "image(preimage(a)) not included in a"


# image(preimage(a)) ⊆ a
preimage_tests = [
    ("ZPre1_invertible",
     ["{(i, j) | 1 <= i <= 50, 1 <= j <= 50}"],
     ["(i, j -> 2i, 3j)"]),

    ("ZPre1_flatten",
     ["{(i, j) | 1 <= i <= 50, 1 <= j <= 50}"],
     ["(i -> i, i)"]),

    ("ZPre1_extend",
     ["{(i, j) | 1 <= i <= 50, 1 <= j <= 50}"],
     ["(i, j, k -> i, j)"]),
]



from do_tests import do_tests
do_tests(preimage_tests, preimage_test_func)
