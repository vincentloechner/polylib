"""
test_ZImage.py — test image and preimage of LBLs

Usage:
    python3 test_ZImage.py
"""
from pypolylib import LBL, Transfo


# test function: check that a in f^{-1}(f(a))
# a ⊆ preimage(image(a))
def image_test_func(a, f):
    img = a.image(f)
    pre = img.preimage(f)
    assert a in pre, "a not included in preimage(image(a))"


# (name, a:lbl, f:transfo)
image_tests = [
    ("ZImage1_preimage_of_image",
     LBL("{(2i, 2i+j) | 1 <= i <= 10, 1 <= j <= 10}"),
     Transfo("(i, j -> 3i+1, 2i+5j)")),

    ("ZImage2_preimage_of_image",
     LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}"),
     Transfo("(i, j -> 2i+1, i+3j)")),

    ("ZImage3_preimage_of_image",
     LBL("{(i, j) | 0 <= i <= 20, 0 <= j <= 20}"),
     Transfo("(i, j -> i+j, i-j)")),

    ("ZImPre1_roundtrip",
     LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}"),
     Transfo("(i, j -> i+1, j+1)")),

    ("ZImPre_flatten",
     LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}"),
     Transfo("(i, j -> i)")),

    ("ZImPre_flatten_pointy",
     LBL("{(i, j) | 0 <= i <= 100, 10i <= 50j <= 15i}"),
     Transfo("(i, j -> i)")),

    ("ZImPre_extend_dim",
     LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}"),
     Transfo("(i, j -> i, j, i+j)")),
]



from do_tests import do_tests
do_tests(image_tests, image_test_func)
