"""
run_tests.py — Run all PyPolyLib tests in python subprocesses

Usage:
    python3 run_tests.py test_*.py
"""
import subprocess, sys

failed = 0
for f in sys.argv[1:]:
    print(f"{'='*15} {f} {'='*15}")
    try:
        r = subprocess.run(["python3", f], timeout=120)
        if r.returncode != 0:
            failed += r.returncode
            print(f"{r.returncode} FAIL in {f}")
    except subprocess.TimeoutExpired:
        print(f"\n\033[31mTIMEOUT in {f}\033[0m")
        failed += 1


print(f"\n{'='*40}")
if failed:
    print(f"\033[31mTotal: {failed} failed\033[0m")
    print(f"{'='*40}")
    sys.exit(1)
print(f"All tests passed")
print(f"{'='*40}")
