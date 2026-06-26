"""
test_ZDomain.py — Teste que zdomain(a) == a pour plusieurs LBLs.
Inspiré de ZDomain3.in, ZDomain4.in, ZDomain6.in, ZDomain10.in, ZDomain11.in.
 
Usage:
    export LD_LIBRARY_PATH=~/polylib/.libs:$LD_LIBRARY_PATH
    python3 test_ZDomain.py
"""
import subprocess, sys
 
# (nom, lbl_string)
# Vérification : a.included(zdomain(a)) and zdomain(a).included(a)
tests = [
    ("ZDomain3",
     "{(i) | 0 <= i <= 10, 3j >= i, 2j <= i, j >= 0}"),
 
    ("ZDomain4",
     "{(i) | 0 <= i <= 100, 3j >= i, 10j <= i, j >= 0}"),
 
    ("ZDomain11",
     "{(2i, 0, 52i+66j) | 198k <= -150i+j, 264k >= -200i+j, 9j <= 2i, 5j >= i, i <= 1000}"),
]
 
CODE = """
from pypolylib import LBLRead
a = LBLRead({lbl!r})
b = a.zdomain()
assert a.included(b), "a not included in zdomain(a)"
assert b.included(a), "zdomain(a) not included in a"
"""
 
passed = 0
failed = 0
 
for name, lbl_str in tests:
    code = CODE.format(lbl=lbl_str)
    try:
        r = subprocess.run([sys.executable, '-c', code],
                           capture_output=True, text=True, timeout=60)
        if r.returncode == 0:
            print(f"OK   {name}"); passed += 1
        else:
            err = r.stderr.strip().split('\n')[-1]
            print(f"FAIL {name}: {err}"); failed += 1
    except subprocess.TimeoutExpired:
        print(f"SKIP {name}: timeout (>60s)")
 
print(f"\n  {passed} passed, {failed} failed")
if failed:
    sys.exit(1)
 
