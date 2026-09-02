# PyPolyLib

This is PyPolyLib, a Python interface to PolyLib.

Written by:
- Vincent LOECHNER <loechner@unistra.fr>
- Siham AIT KHORSA <siham.ait-khorsa@etu.unistra.fr> and
  Sabrina BOUGARECH <sabrina.bougarech@etu.unistra.fr>
  (during their 2026 internship under the supervision of Vincent LOECHNER)

## Requirements

- PolyLib (included) with GMP support
- CMake
- Python3
- pybind11
- NumPy, SciPy and PyVista (optional, for plotting)

PyPolyLib requires the `pybind11` module to build the C++ Python bindings
to the library.

The `numpy`, `scipy` and `pyvista` modules are used by the plotting
functionality, but they are optional (an `ImportError` will be raised if you
try to plot an LBL when they are not installed).

## Installation
0. Get PolyLib from the git:
```sh
git clone https://github.com/vincentloechner/polylib
cd polylib
./autogen.sh
```

1. Build PolyLib

Build the GMP version of PolyLib from the main PolyLib directory in
a `build/` subdirectory with:
```sh
# From the PolyLib main directory:
mkdir build && cd build

# on recent MacOS you need to configure like this
# (to find gmp and configure the generated uninstalled libs rpath correctly)
# on Linux with gmp installed system-wide no configure options are needed.
# add the --prefix option to choose your installation directory (optional,
# default installation on /usr/local).
../configure \
  --with-libgmp=$(pkg-config gmp --variable=prefix) \
  --libdir=$(pwd)/.libs

make -j 20              # build PolyLib

make -j 20 check        # optional

# For testing you can directly use the PolyLib build, you do not need to
# install it system-wide.
# If you want to install it on a longer-term basis run:
make -j 20 install      # optional, use sudo if necessary

# return to the PolyLib main directory
cd ..
```

2. Enter PyPolyLib and prepare the Python environment: set a Python venv
(optional) and install all required Python packages:
```sh
cd PyPolyLib
python3 -m venv ./polylib_venv && source ./polylib_venv/bin/activate
python3 -m pip install pybind11 numpy pyvista scipy
```

3. Build pypolylib:

- using the local build polylib located in `../build` with:
```sh
cmake -S . -B build/polylib-local \
  -DPOLYLIB_ROOT=../build \
  -DPOLYLIB_USE_UNINSTALLED=ON
cmake --build build/polylib-local
```

- or using the installed polylib (if you installed it system-wide) with:
```sh
cmake -S . -B build/polylib-system
cmake --build build/polylib-system
```
If it is not found in a standard path, you can specify where to find the
pkg-config `polylibgmp.pc` file using the PKG_CONFIG_PATH shell variable.

To test the PyPolyLib python build, you just need to set your
`PYTHONPATH` shell variable:
```sh
export PYTHONPATH=$(pwd)/build/polylib-local
# or:
export PYTHONPATH=$(pwd)/build/polylib-system
```


## Usage

Do not forget to set your environment variable `PYTHONPATH` and activate
your python venv (located in the PyPolyLib/ directory if you followed the
above instructions), especially if you run another shell.

### Basic syntax
- Import the library using:
```py
from pypolylib import LBL, Transfo
```

- Read an LBL from a string, for example:
```py
a = LBL("{(2i, 2i+j) | 1 <= i <= 5, 1 <= j <= i}")
```
The left-hand side expression is a vector of linear functions of any number
of variables. The right-hand side is a combination of linear inequalities and
equalities, any supplementary variable appearing here will be considered as
an existential variable.
The braces and the parenthesis are optional.

- Available operations in a nutshell:
```py
a + b == a.union(b)
a - b == a.difference(b)
a * b == a.intersection(b)
a in b == a.included(b)     # geometric LBL inclusion a ⊆ b
(1,2,3) in a == a.contains_point((1,2,3)) # "pythonic" inclusion
for z in a:
  ...  # scan all points z in a
set(a) # builds the python set of points in a
(a == b) == ((a in b) and (b in a))
a.zdomain()  # computes the corresponding Z-domain, eliminating all
             # existential variables (can be complex!)
a.disjoint() # computes a disjoint union representation of
             # the union of LBLs in a (can be complex!)
a.plot() # plots the LBL in a window (see help for available options)
a.sLBL_list() # get a Python list of single LBLs from the union of LBLs
              # represented by a
```

- Transformation operations:
```py
f = Transfo("(i,j -> 3i+2j+1)") # define a transformation function
g = Transfo("(i,j -> j,i)")
f * g == f.compose(g)
f.inverse()   # *integer* inverse of this function (see help)
f(a) == a.image(f) # image of the LBL a by f
a.preimage(f) # computes the preimage of LBL a by function f.
              # Notice that if the preimage by f is not consistent to form
              # an LBL only its LBL part is returned, which means that:
              # f(a.preimage(f)) is not necessarily equal to a (but it is
              # included in a).
```


### Python program examples
The `./examples/` directory contains some example python programs.

### Interactive example
You can get a simple example running in an interactive python
console by typing:
`python3 -i examples/demo.py`

It runs the following demo (with some extra pretty printing):
```py
>>> from pypolylib import LBL
>>> a = LBL("{(2i, 2i+j) | 1 <= i <= 5, 1 <= j <= i}")
>>> b = LBL("{(i,j) | 0 <= i <= 6, 0 <= j <= 8}")
>>> print(f"a = {a}")
a = LBL("{(2i, j) | i <= 5, -2i+j - 1 >= 0, 3i-j >= 0}")
>>> print(f"b = {b}")
b = LBL("{(i, j) | i >= 0, i <= 6, j >= 0, j <= 8}")
>>>
>>> print(f"a + b = {a + b}")
a + b = LBL("{(2i, j) | i <= 5, -2i+j - 1 >= 0, 3i-j >= 0}") + \
LBL("{(i, j) | i >= 0, i <= 6, j >= 0, j <= 8}")
>>> print(f"a - b = {a - b}")
a - b = LBL("{(2i, j) | i <= 5, -2i+j - 1 >= 0, 3i-j >= 0, j >= 9}")
>>> print(f"a * b = {a * b}")
a * b = LBL("{(2i, j) | -2i+j - 1 >= 0, 3i-j >= 0, i <= 3, j <= 8}")
>>>
>>> for z in a:
...   print(z)
...
(2, 3)
(4, 5)
(4, 6)
(6, 7)
(6, 8)
(6, 9)
(8, 9)
(8, 10)
(8, 11)
(8, 12)
(10, 11)
(10, 12)
(10, 13)
(10, 14)
(10, 15)
>>>
>>> print(f"set(a) = {set(a)}")
set(a) = {(10, 15), (10, 11), (10, 14), (6, 8), (4, 6), (8, 10), (10, 13), (2, 3), (6, 7), (4, 5), (8, 9), (10, 12), (8, 12), (6, 9), (8, 11)}
>>>
>>> # plot a:
>>> a.plot()
>>>
>>> # nice plot of a, a-b, b-a and i=a*b in two separate windows:
>>> i = a * b
>>> # plot (the border of a) and (a as the two subsets: (a-i) + i) in window1:
>>> a.plot(subplot=True, show_points=False, color="black")
>>> (a-i).plot(subplot=True, color="blue")
>>> i.plot(interactive_update=True, color="cyan")
>>> # plot (the border of b) and (b-i) in window2:
>>> b.plot(subplot=True, show_points=False, color="black")
>>> (b-i).plot(color="green")
>>>
```
