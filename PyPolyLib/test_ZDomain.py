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
    ("ZDomain2b",
     "{(i) | 0 <= i, 30j >= i, 29j <= i}"),

    ("ZDomain3 ",
     "{(i) | 0 <= i <= 10, 3j >= i, 2j <= i}"),

    ("ZDomain4 ",
     "{(i) | 0 <= i <= 100, 3j >= i, 10j <= 5i-1, 11i <= 20j + 100}"),

    ("ZDomain6 ",
     "{(2i, 8i+11j) | 24i+44j-k+40>=0, -18i-33j+k-30>=0, -18i-33j+k-30>=0, 24i+44j-l+40>=0, -i+5l>=0, 2i-9l>=0, 2i-9k>=0, -i+5k>=0, i>=51}"),

    ("ZDomain10",
     "{(i, j) | i>=101,  3i-3j+k>=0, -4i+4j-k>=0, -i+10k>=0, i-9k>=0}"),

    ("ZDomain12 (the Omega test nightmare)",
     "{(i, j) | 11i+13j>=27, 11i+13j<=45, 7i-9j+10>=0, -7i+9j+4>=0}"),
 
    # very long one:
    # ("ZDomain11",
    #  "{(2i, 0, 52i+66j) | 198k <= -150i+j, 264k >= -200i+j, 9j <= 2i, 5j >= i, i <= 1000}"),
]
 
from pypolylib import LBL
 
passed = 0
failed = 0

for name, lbl_str in tests:
    try:
        print(f"{name}... ", end="", flush=True)
        a = LBL(lbl_str)
        b = a.zdomain()
        # assert a.included(b), "a not included in zdomain(a)"
        # assert b.included(a), "zdomain(a) not included in a"
        assert a == b
        print("OK")
        passed += 1
    except AssertionError:
        print("FAIL")
        failed += 1
 
print(f"\n  {passed} passed, {failed} failed")
if failed:
    sys.exit(1)
