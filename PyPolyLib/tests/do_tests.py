import sys
import traceback
from pypolylib import LBL, Transfo, PolylibClose

def do_nothing(*args):
    pass

def run(tests, test_func, fail_func):

    passed = 0
    failed = 0
    for name, *test_args in tests:
        print(f"{name:22s} ", end="", flush=True)
        try:
            # init LBLs
            lbl = list(LBL(s) for s in test_args[0])
            trans = []

            if len(test_args) == 2:
                # init tranformation matrices
                trans = list(Transfo(s) for s in test_args[1])

        except:
            print("\033[31mFAIL\033[0m: cannot read LBL(s)")
            print(traceback.format_exc())
            failed += 1
            continue

        try:
            test_func(*lbl, *trans)
            print("OK")
            passed += 1
        except:
            print("\033[31mFAIL\033[0m")
            print(traceback.format_exc())
            fail_func(*lbl)
            failed += 1

    print(f"{passed} passed, {failed} failed")
    return failed


def do_tests(tests, test_func, fail_func=do_nothing):
    failed = run(tests, test_func, fail_func)

    # clean buffers
    PolylibClose()

    # exit code = number of failed tests (0 == success)
    sys.exit(failed)
