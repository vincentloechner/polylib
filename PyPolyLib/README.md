# PyPolyLib

This is PyPolyLib, a Python interface to PolyLib.

Written by:
- Vincent LOECHNER <loechner@unistra.fr>
- Siham AIT KHORSA <siham.ait-khorsa@etu.unistra.fr> and
  Sabrina BOUGARECH <sabrina.bougarech@etu.unistra.fr>
  (under the supervision of Vincent LOECHNER <loechner@unistra.fr>)


## Installation

First build the GMP version of PolyLib from the main directory in a
`build/` subdirectory (optional):
```sh
# From the PolyLib main directory:
mkdir build && cd build
# on MacOS you need something like this (no options needed on Linux with gmp installed):
../configure --with-libgmp=$(pkg-config gmp --variable=prefix) --libdir=$(pwd)/.libs
make -j 20
make -j 20 check
# You can directly use the polylib build, you do not need to install it system-wide.
# Run `make install` if you want to install it.

# Enter the PyPolyLib directory
cd ../PyPolyLib
```

Prepare the Python installation: set a Python venv (optional), get all
required Python packages, and set `PYTHONPATH`.
```sh
cd ../PyPolyLib
python3 -m venv polylib && source ./polylib/bin/activate
python3 -m pip install pybind11 numpy matplotlib shapely
export PYTHONPATH="$(pwd):$PYTHONPATH"
```

Build pypolylib_core from the current directory using:
```sh
make POLYLIB_BUILD=../build
````


## Usage

Do not forget to set your environment variables `PYTHONPATH` and
`LD_LIBRARY_PATH`, and activate your python venv.

### Python program examples
Run the examples in the `./tests/` and `./examples/` directories.

### Interactive example
Run `python3` and type:
```py
>>> from pypolylib import LBL
>>> a = LBL("{(2i, 2i+j) | 1 <= i <= 5, 1 <= j <= i}")
>>> b = LBL("{(i,j)|j<5}")
>>> a
{(2i, j) | i <= 5, -2i+j - 1 >= 0, 3i-j >= 0}
>>> b
{(i, j) | j <= 4, 0 + 1 >= 0}
>>> # union
>>> a + b
{(2i, j) | i <= 5, -2i+j - 1 >= 0, 3i-j >= 0} UNION
{(i, j) | j <= 4, 0 + 1 >= 0}
>>> # difference
>>> a - b
{(2i, j) | i <= 5, -2i+j - 1 >= 0, 3i-j >= 0, j >= 5}
>>> # intersection
>>> a * b
{(2i, j) | -2i+j - 1 >= 0, 3i-j >= 0, j <= 4}
>>> 
>>> # scan (bounded LBL)
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
>>> set(a)
{(10, 15), (10, 11), (10, 14), (6, 8), (4, 6), (8, 10), (10, 13), (2, 3), (6, 7), (4, 5), (8, 9), (10, 12), (8, 12), (6, 9), (8, 11)}
>>> 
>>> # visualization window (bounded 2D or 3D LBL)
>>> a.plot()
>>> b = LBL("{(i,j)|0<=i<=5,0<=j<=5}")
>>> (a+b).plot()
>>> 
```
