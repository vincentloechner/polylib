# quick and dirty makefile, just to test

POLYLIB_BUILD=../../../build/
PYB_INC=$(shell python3 -m pybind11 --includes)
GMP_INC=$(shell pkg-config --cflags gmp)
GMP_LIB=$(shell pkg-config --libs gmp)
PY_LIB=$(shell python3-config --embed --ldflags)
PY_EXT=$(shell python3 -m pybind11 --extension-suffix)


# $ c++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup $(python3 -m pybind11 --includes) example.cpp -o example$(python3-config --extension-suffix)

all: pypolylib_core$(PY_EXT)


pypolylib_core$(PY_EXT): pypolylib_core.cpp
# this one works on MacOS:
	g++ -O3 -Wall -shared -std=c++11  -undefined dynamic_lookup -fPIC -I../../../include ${GMP_INC} ${PYB_INC} $< -o $@ -L${POLYLIB_BUILD}/.libs -lpolylibgmp ${GMP_LIB} ${PY_LIB}

# this one is ok on Linux with everything installed correctly
# c++ -Wall -shared -std=c++11 -fPIC $< ../../../.libs/libpolylibgmp.so $(PYB_INC) -I../../../include -o $@
