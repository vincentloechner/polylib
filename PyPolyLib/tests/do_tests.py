import sys
import traceback
from pypolylib import LBL, Transfo

def do_nothing(*args):
    pass

def do_tests(tests, test_func, fail_func=do_nothing):
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

    # exit code = number of failed tests (0 == success)
    sys.exit(failed)
