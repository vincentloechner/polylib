# quick calls to cmake/python to build & install
#
# requires Python's build module:
# python3 -m pip install --break-system-packages build
# (option --break-system-packages can be avoided using a personnal venv)

.PHONY: all local install check clean

CMAKE=cmake
PYTHON=python3

# build the python module using the local polylib build
all: build/polylib-local
	$(CMAKE) --build build/polylib-local
	$(PYTHON) -m build --wheel  -Ccmake.define.POLYLIB_USE_UNINSTALLED=ON

build/polylib-local:
	# initial cmake to create the build tree
	$(CMAKE) -S . -B build/polylib-local -DPOLYLIB_DIR=../build  -DPOLYLIB_USE_UNINSTALLED=ON

install: all
	$(PYTHON) -m venv test-venv
	# just update the pypolylib module, not all the standard ones:
	test-venv/bin/pip uninstall -y pypolylib
	test-venv/bin/pip install dist/*.whl

check:
	test-venv/bin/python tests/run_all_tests.py


# # this was a quick and dirty makefile, just to test
# # (don't forget to set LD_LIBRARY_PATH to find the lib at runtime)
# POLYLIB_BUILD?=../build/
# POLYLIB=$(shell pkg-config -cflags -libs ${POLYLIB_BUILD}/polylibgmp-uninstalled.pc || pkg-config -cflags -libs polylibgmp || echo "polylib_NOTFOUND")
# PYBIND= $(shell python3 -m pybind11 --includes)
# PY_LIB= $(shell python3-config --embed --ldflags)
# PY_EXT= $(shell python3-config --extension-suffix)
# # this is required if gmp is not installed in a standard path:
# GMP=    $(shell pkg-config --cflags -libs-only-L gmp)
#
# all: pypolylib_core${PY_EXT}
#
# pypolylib_core${PY_EXT}: pypolylib_core.cpp
# 	g++ -O3 -Wall -shared -std=c++11 -fPIC -I../include $< -o $@ ${POLYLIB} ${GMP} ${PYBIND} ${PY_LIB}
# clean:
# 	rm pypolylib_core${PY_EXT}
# check:
# 	python3 tests/run_all_tests.py
