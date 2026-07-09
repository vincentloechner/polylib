# quick and dirty makefile, just to test

# trying built uninstalled version of polylib
# (don't forget to set LD_LIBRARY_PATH to find the lib at runtime)
POLYLIB_BUILD?=../build/
POLYLIB=$(shell pkg-config -cflags -libs ${POLYLIB_BUILD}/polylibgmp-uninstalled.pc || pkg-config -cflags -libs polylibgmp || echo "polylib_NOTFOUND")
PYBIND= $(shell python3 -m pybind11 --includes)
PY_LIB= $(shell python3-config --embed --ldflags)
PY_EXT= $(shell python3-config --extension-suffix)
# this is required if gmp is not installed in a standard path:
GMP=    $(shell pkg-config --cflags -libs-only-L gmp)


all: pypolylib_core${PY_EXT}

pypolylib_core${PY_EXT}: pypolylib_core.cpp
	g++ -O3 -Wall -shared -std=c++11 -fPIC -I../include $< -o $@ ${POLYLIB} ${GMP} ${PYBIND} ${PY_LIB}

clean:
	rm pypolylib_core${PY_EXT}

check:
	python3 tests/run_tests.py tests/test_*
