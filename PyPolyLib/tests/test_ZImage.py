"""
test_ZImage.py — test image and preimage of LBLs

Usage:
    python3 test_ZImage.py
"""
from pypolylib import LBL, Transfo


# test function: check that f^{-1}(f(a)) \in a
# preimage(image(a)) ⊆ a
def image_test_func(a, f):
    img = a.image(f)
    pre = img.preimage(f)
    assert a in pre, "a not included in preimage(image(a))"


# (name, [lbl_str, ...], [transfo_str, ...])
image_tests = [
    ("ZImage1_preimage_of_image",
     ["{(2i, 2i+j) | 1 <= i <= 10, 1 <= j <= 10}"],
     ["(i, j -> 3i+1, 2i+5j)"]),

    ("ZImage2_preimage_of_image",
     ["{(i, j) | 1 <= i <= 10, 1 <= j <= 10}"],
     ["(i, j -> 2i+1, i+3j)"]),

    ("ZImage3_preimage_of_image",
     ["{(i, j) | 0 <= i <= 20, 0 <= j <= 20}"],
     ["(i, j -> i+j, i-j)"]),

    ("ZImPre1_roundtrip",
     ["{(i, j) | 1 <= i <= 10, 1 <= j <= 10}"],
     ["(i, j -> i+1, j+1)"]),

    ("ZImPre_flatten",
     ["{(i, j) | 1 <= i <= 10, 1 <= j <= 10}"],
     ["(i, j -> i)"]),

    ("ZImPre_flatten_pointy",
     ["{(i, j) | 0 <= i <= 100, 10i <= 50j <= 15i}"],
     ["(i, j -> i)"]),

    ("ZImPre_extend_dim",
     ["{(i, j) | 1 <= i <= 10, 1 <= j <= 10}"],
     ["(i, j -> i, j, i+j)"]),
]



from do_tests import do_tests
do_tests(image_tests, image_test_func)
