# PyPolyLib

This is PyPolyLib, a Python interface to PolyLib.

Written by:
- Siham AIT KHORSA <siham.ait-khorsa@etu.unistra.fr>
- Sabrina BOUGARECH <sabrina.bougarech@etu.unistra.fr>

under the supervision of:
- Vincent LOECHNER <loechner@unistra.fr>


## Installation

Install the required Python packages:
```sh
# optional: create your python venv with
python3 -m venv pplenv && source ./pplenv/bin/activate
python3 -m pip install pybind11 numpy matplotlib shapely
```

Build pypolylib_core from the current directory using:
```sh
make
````
You may need to adjust some environment variables to map your OS and PolyLib installation, see the first lines of the Makefile.

## Usage

### Interactive example
Run `python3`in the current directory and then:
```py
>>> from pypolylib import LBL
>>> a = LBL("{(2i, 2i+j) | 1 <= i <= 5, 1 <= j <= i}")
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
>>> # open a visualization window
>>> a.plot()
>>> 
>>> # OPERATIONS:
>>> b = LBL("{(i,j)|j<5}")
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
```
