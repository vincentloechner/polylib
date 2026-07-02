"""
test_ZImage.py — Teste image et préimage de LBLs sous des transformations.
Inspiré de ZImage1.in ... ZImage4.in, ZImPre1.in, ZPre1.in.

Usage:
    export LD_LIBRARY_PATH=~/polylib/.libs:$LD_LIBRARY_PATH
    python3 test_ZImage.py
"""
import subprocess, sys

# (nom, lbl_a, transfo, type)
# type "image"   : preimage(image(a)) ⊆ a
# type "preimage": image(preimage(a)) ⊆ a
tests = [
    ("ZImage1_preimage_of_image",
     "{(2i, 2i+j) | 1 <= i <= 10, 1 <= j <= 10}",
     "(i, j -> 3i+1, 2i+5j)",
     "image"),

    ("ZImage2_preimage_of_image",
     "{(i, j) | 1 <= i <= 10, 1 <= j <= 10}",
     "(i, j -> 2i+1, i+3j)",
     "image"),

    ("ZImage3_preimage_of_image",
     "{(i, j) | 0 <= i <= 20, 0 <= j <= 20}",
     "(i, j -> i+j, i-j)",
     "image"),

    ("ZPre1_image_of_preimage",
     "{(i, j) | 1 <= i <= 50, 1 <= j <= 50}",
     "(i, j -> 2i, 3j)",
     "preimage"),

    ("ZImPre1_roundtrip",
     "{(i, j) | 1 <= i <= 10, 1 <= j <= 10}",
     "(i, j -> i+1, j+1)",
     "image"),
]

CODE_IMAGE = """
from pypolylib import LBLRead, Transfo
a = LBLRead({a!r})
f = Transfo({f!r})
img = a.image(f)
pre = img.preimage(f)
assert pre.included(a), "preimage(image(a)) not included in a"
"""

CODE_PREIMAGE = """
from pypolylib import LBLRead, Transfo
a = LBLRead({a!r})
f = Transfo({f!r})
pre = a.preimage(f)
img = pre.image(f)
assert img.included(a), "image(preimage(a)) not included in a"
"""

passed = 0
failed = 0

for name, lbl_a, transfo, typ in tests:
    if typ == "image":
        code = CODE_IMAGE.format(a=lbl_a, f=transfo)
    else:
        code = CODE_PREIMAGE.format(a=lbl_a, f=transfo)
    r = subprocess.run([sys.executable, '-c', code],
                       capture_output=True, text=True, timeout=30)
    if r.returncode == 0:
        print(f"OK   {name}"); passed += 1
    else:
        err = r.stderr.strip().split('\n')[-1]
        print(f"FAIL {name}: {err}"); failed += 1

print(f"\n  {passed} passed, {failed} failed")
if failed:
    sys.exit(1)