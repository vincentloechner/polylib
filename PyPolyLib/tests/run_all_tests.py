"""
run_tests.py — Run all PyPolyLib tests in python subprocesses

Run the test_*.py files from the directory where run_tests.py lies.
Timeout can be set with the "-t SECONDS" option

Usage:
    python3 run_tests.py
"""
import subprocess, sys
from pathlib import Path
import argparse

# parse timeout arg:
parser = argparse.ArgumentParser(description="Test runner")
parser.add_argument(
    "-t", "--timeout",
    type=int,
    default=60,
    help="Timeout (seconds) per test (default: 60)"
)
args = parser.parse_args()
timeout = args.timeout

# get file names
script_dir = Path(__file__).parent
files = sorted(script_dir.glob("test_*.py"), key=lambda f: f.name.lower())

# run tests
failed = 0
for f in files:
    print(f"{'='*15} {f} {'='*15}")
    try:
        r = subprocess.run(["python3", f], timeout=timeout)
        if r.returncode != 0:
            failed += r.returncode
            print(f"{r.returncode} FAIL in {f}")
    except subprocess.TimeoutExpired:
        print(f"\n\033[31mTIMEOUT in {f}\033[0m")
        failed += 1

# final message
print(f"\n{'='*40}")
if failed:
    print(f"\033[31mTotal: {failed} failed\033[0m")
    print(f"{'='*40}")
    sys.exit(1)
print(f"All tests passed")
print(f"{'='*40}")
