# quick and dirty makefile, just to test
PYB_INC=$(shell python3 -m pybind11 --includes)
PY_EXT=$(shell python3-config --extension-suffix)

all: pypolylib_core$(PY_EXT)

pypolylib_core$(PY_EXT): pypolylib_core.cpp
	c++ -Wall -shared -std=c++11 -fPIC $< ../../../.libs/libpolylibgmp.so $(PYB_INC) -I../../../include -o $@
