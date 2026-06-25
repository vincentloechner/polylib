"""
run_tests.py — Lance tous les tests PyPolyLib dans des sous-processus séparés.

Usage:
    export LD_LIBRARY_PATH=~/polylib/.libs:$LD_LIBRARY_PATH
    python3 run_tests.py
"""
import subprocess, sys

tests = {
    "ZAlloc3d":
        'from pypolylib import LBLRead; a=LBLRead("{(i+j)|0<=i<=10,i=2j}"); e=LBLRead("{(3i)|0<=i<=5}"); assert a==e',
    "ZDiffInter6a":
        'from pypolylib import LBLRead; a=LBLRead("{(i,j)|1<=i<=20,1<=j<=20}"); b=LBLRead("{(2i,j)|1<=i<=50,j=10}"); assert (a-b+a*b)==a',
    "ZDomain3":
        'from pypolylib import LBLRead; a=LBLRead("{(i)|0<=i<=10,3j>=i,2j<=i,j>=0}"); b=a.zdomain(); assert a.included(b) and b.included(a)',
    "ZImage1":
        'from pypolylib import LBLRead,Transfo; a=LBLRead("{(2i,2i+j)|1<=i<=10,1<=j<=10}"); f=Transfo("(i,j->3i+1,2i+5j)"); assert a.image(f).preimage(f).included(a)',
    "ZUnion":
        'from pypolylib import LBLRead; a=LBLRead("{(i,j)|1<=i<=5,1<=j<=5}"); b=LBLRead("{(i,j)|8<=i<=10,8<=j<=10}"); u=a+b; assert a.included(u) and b.included(u)',
    "incompatible_dimensions":
        'from pypolylib import LBLRead\na=LBLRead("{(i,j)|0<=i<=5,0<=j<=i}")\nc=LBLRead("{(i,j,k)|0<=i<=5,0<=j<=i,0<=k<=i}")\nfor op in [lambda:a+c,lambda:a*c,lambda:a-c]:\n  try: op(); assert False\n  except ValueError: pass',
    "unbounded":
        'from pypolylib import LBLRead\nc=LBLRead("{(i,j)|0<=i<=5,j>=0}")\ntry: list(c); assert False\nexcept ValueError: pass',
}

passed = 0
failed = 0
for name, code in tests.items():
    r = subprocess.run([sys.executable, '-c', code],
                       capture_output=True, text=True, timeout=30)
    if r.returncode == 0:
        print(f"OK   {name}"); passed += 1
    else:
        err = r.stderr.strip().split('\n')[-1]
        print(f"FAIL {name}: {err}"); failed += 1

print(f"\n{'='*40}")
print(f"  {passed} passed, {failed} failed")
print(f"{'='*40}")
if failed: sys.exit(1)
