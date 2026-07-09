"""
test_ZDomain.py — Check if zdomain(a) == a for several example LBLs.
"""
from pypolylib import LBL


# check if: a.included(zdomain(a)) and zdomain(a).included(a)
def test_func(a):
    b = a.zdomain()
    assert a.included(b), "a not included in zdomain(a)"
    assert b.included(a), "zdomain(a) not included in a"
    # equivalent to a == b, with more precise output

def fail_func(a):
    print("a = {a}")
    print(f"a.zdomain() = {a.zdomain()}")

# list of (name:str, [a:lbl])
tests = [
    ("empty",
     ["{(i) | 1<0}"]),

    ("ZDomain2b",
     ["{(i) | 0 <= i, 30j >= i, 29j <= i}"]),

    ("ZDomain3",
     ["{(i) | 0 <= i <= 10, 3j >= i, 2j <= i}"]),

    ("ZDomain4",
     ["{(i) | 0 <= i <= 100, 3j >= i, 10j <= 5i-1, 11i <= 20j + 100}"]),

    ("ZDomain6",
     ["{(2i, 8i+11j) | 24i+44j-k+40>=0, -18i-33j+k-30>=0, -18i-33j+k-30>=0, 24i+44j-l+40>=0, -i+5l>=0, 2i-9l>=0, 2i-9k>=0, -i+5k>=0, i>=51}"]),

    ("ZDomain10",
     ["{(i, j) | i>=101,  3i-3j+k>=0, -4i+4j-k>=0, -i+10k>=0, i-9k>=0}"]),

    # from Pugh91 (the Omega test nightmare)
    ("ZDomain12",
     ["{(i, j) | 11i+13j>=27, 11i+13j<=45, 7i-9j+10>=0, -7i+9j+4>=0}"]),
]

# pretty long ones below, not tested by default:
longtests = [
    # this one is from Rui-Juan Jing, Marc Moreno Maza, et al.,
    # Quantifier Elimination Over the Integers
    # https://dl.acm.org/doi/pdf/10.1145/3747199.3747580
    ("ZDomain13",
     ["{(x1, x2) | -98877x1 - 189663x2 - 1798x3 <= 705915, -10109x1 - 5958x2 - 14601x3 <= 31333, -5405x1 + 4965x2 + 3870x3 <= 4303504, 729x1 - 117x2 + 350x3 <= 4561, 677x1 + 465x2 - 540x3 <= 3489}"]),

    ("ZDomain11",
     ["{(2i, 0, 52i+66j) | 198k <= -150i+j, 264k >= -200i+j, 9j <= 2i, 5j >= i, i <= 1000}"]),
]


from do_tests import do_tests
do_tests(tests, test_func, fail_func)

# Can run those longer ones (comment out the previous one)
# but run_tests.py sets a 2 minutes timeout, likely to expire
# do_tests(tests + longtests, test_func, fail_func)
