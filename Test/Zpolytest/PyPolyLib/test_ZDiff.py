"""
test_ZDiff.py — Teste que (a-b) + (a∩b) == a pour plusieurs paires de LBLs.
Inspiré de ZDiff1.in ... ZDiff13.in et ZDiffInter6a.in.

Usage:
    export LD_LIBRARY_PATH=~/polylib/.libs:$LD_LIBRARY_PATH
    python3 test_ZDiff.py
"""
import subprocess, sys

# (nom, lbl_a, lbl_b)
# Vérification : (a-b) + (a∩b) == a
tests = [
    ("ZDiff1",
     "{(i, j) | 1 <= i <= 10, 1 <= j <= 10}",
     "{(i, j) | 5 <= i <= 15, 5 <= j <= 15}"),

    ("ZDiff2",
     "{(i, j) | 0 <= i <= 5, 0 <= j <= i}",
     "{(i, j) | 5 <= i <= 7, 0 <= j <= i}"),

    ("ZDiff3",
     "{(i, j) | 1 <= i <= 20, 1 <= j <= 20}",
     "{(i, j) | 10 <= i <= 30, 10 <= j <= 30}"),

    ("ZDiffInter6a",
     "{(i, j)  | 1 <= i <= 20, 1 <= j <= 20}",
     "{(2i, j) | 1 <= i <= 50, j = 10}"),

    ("ZDiff_triangle",
     "{(i, j) | 0 <= i <= 8, 0 <= j <= i}",
     "{(2i, j) | 0 <= i <= 4, 0 <= j <= 2i}"),
]

CODE = """
from pypolylib import LBLRead
a = LBLRead({a!r})
b = LBLRead({b!r})
inter = a * b
diff  = a - b
u     = diff + inter
assert u == a, f"(a-b)+(a∩b) != a"
"""

passed = 0
failed = 0

for name, lbl_a, lbl_b in tests:
    code = CODE.format(a=lbl_a, b=lbl_b)
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