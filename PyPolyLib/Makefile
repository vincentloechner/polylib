# quick and dirty makefile, just to test

POLYLIB_BUILD=../build/
POLYLIB=$(shell pkg-config -cflags -libs --with-path=${POLYLIB_BUILD} polylibgmp)
GMP=$(shell pkg-config --cflags -libs gmp)
PYBIND=$(shell python3 -m pybind11 --includes)
PY_LIB=$(shell python3-config --embed --ldflags)
PY_EXT=$(shell python3-config --extension-suffix)

# not necessary:
# MACOS=$(shell if [ `uname` = "Darwin" ]; then echo "-undefined dynamic_lookup"; fi)

all: pypolylib_core$(PY_EXT)


pypolylib_core$(PY_EXT): pypolylib_core.cpp
	g++ -O3 -Wall -shared -std=c++11 ${MACOS} -fPIC -I../include $< -o $@ ${POLYLIB} ${GMP} ${PYBIND} ${PY_LIB}
