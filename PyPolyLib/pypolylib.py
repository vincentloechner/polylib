"""
pypolylib.py — Python interface to the PolyLib library (LBLs).

This module exposes the LBL (image of Z-polyhedron) and Transfo (matrix)
types for manipulating sets of integer points in Python.

Typical usage:

```
from pypolylib import LBL, Transfo

# Create LBLs
a = LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
b = LBL("{(2i, j) | 1 <= i <= 50, 5 <= j, j < 10}")

#  Set operations
inter = a * b      # intersection
diff  = a - b      # difference
union = a + b      # union

#  Check for inclusion
print(a.included(inter + diff)) # True or False

#  Check for equality (mutual inclusion)
print(a == b)                   # True or False

# Image under a transformation
f = Transfo("(i,j -> 2i+1, i+3j)")
print(a.image(f))               # image of a under f
print(a.preimage(f))            # preimage of a under f

# Composing and inverting transformations
g = Transfo("(i,j -> i+1, 2j)")
h = f * g                       # composition
fi = f.inverse()                # integer inverse

# Iterate over LBL points (bounded LBL)
for pt in a:
  print(pt)                     # print each point (as a tuple of integers)

# get the set of points (bounded LBL, as a python set)
s = set(a)

# Plot (bounded LBL only)
a.plot()
```
"""


import math
import gc

import pypolylib_core as pl
import lbl_io as io
from lbl_plot import lbl_plot


# ──────────────────────────────────────────────────────────────
# Functions
# ──────────────────────────────────────────────────────────────
def LBLRead(s):
  """
  Creates an LBL from a symbolic string.

  This is a shortcut for LBL(s).

  Args:
    s (str): String in the format “{(expr1, ...) | constraints}”

  Returns:
    LBL: The corresponding LBL object

  Exemple:
    >>> a = LBLRead("{(i+j) | 0 <= i <= 10, i = 2j}")
    >>> print(a)
    {(3i) | i >= 0, i <= 5}
  """
  return LBL(s)


def PolylibClose():
  """Cleanup memory (PolyLib allocates reusable buffers).
  """
  gc.collect()
  pl.PolylibClose()


# ──────────────────────────────────────────────────────────────
# LBL Class
# ──────────────────────────────────────────────────────────────
class LBL:
  """
  An LBL is a union of affine image of Z-domains — set of integer points.

  It is stored as a pybind11 polylib _lbl:
  pl.LBL: {Matrix* Lat, Polyhedron* P, pl.LBL* next}

  Creation :
    a = LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")

  Operations :
    a + b   # union
    a - b   # difference
    a * b   # intersection
    etc.
  """
  _lbl:pl.LBL = None
  def __init__(self, obj=None):
    """
    Initializes an LBL, optionally from a symbolic string.

    Arguments:
      obj (str, optional): string in the format “{(expr) | constraints}”.
      obj (pl.LBL, optional): pl.LBL object.
      If None, creates an empty LBL.
    """
    if obj is not None:
      if isinstance(obj, str):
        self._lbl = io.LBLRead(obj)
      elif isinstance(obj, pl.LBL):
        self._lbl = obj
      else:
        raise TypeError(f"expected a string or a pypolylib_core.LBL object, "
                        f" got <{type(obj).__name__}>")

  def __repr__(self):
    """
    Returns a string representation of the LBL, e.g., {(3i) | i >= 0, i <= 5}
    """
    if self._lbl is None:
      return "None"
    return io.LBLRepr(self._lbl)

  ######################## operations ########################
  def _check_compat(self, other):
    if not isinstance(other, LBL):
      raise TypeError(f"expected an LBL, got <{type(other).__name__}>")
    if self._lbl.Lat.nbrows != other._lbl.Lat.nbrows:
      raise ValueError(f"incompatible dimensions: {self._lbl.Lat.nbrows - 1} "
                       f"vs {other._lbl.Lat.nbrows - 1}")

  def intersection(self, other):
    """
    Calculates the intersection of this LBL with another one.

    Args:
      other (LBL): The other LBL (same dimension)

    Returns:
      LBL: The intersection of the two LBLs
    """
    self._check_compat(other)
    return LBL(self._lbl.intersection(other._lbl))

  def difference(self, other):
    """
    Calculates the difference between this LBL and another one (self - other).

    Args:
      other (LBL): The other LBL (same dimension)

    Returns:
      LBL: The difference between the two LBLs
    """
    self._check_compat(other)
    return LBL(self._lbl.difference(other._lbl))

  def union(self, other):
    """
    Calculates the union of this LBL with another one.

    Args:
      other (LBL): The other LBL (same dimension)

    Returns:
      LBL: The union of the two LBLs
    """
    self._check_compat(other)
    return LBL(self._lbl.union(other._lbl))

  def included(self, other):
    """
    Checks whether this LBL is contained within another one.

    Args:
      other (LBL): The other LBL

    Returns:
      bool: True if self ⊆ other, False otherwise
    """
    self._check_compat(other)
    return self._lbl.included(other._lbl)

  def contains_point(self, point):
    """
    Checks whether this LBL contains a point (iterable of correct dimension)

    Args:
      point (iterable)

    Returns:
      bool: True if the lbl contains point, False otherwise
    """
    # check if point is an iterable (not an str or bytes) or raise a type error
    iter(point)
    if isinstance(point, (str, bytes)):
      raise TypeError(f"expected a vector of values, "
                      f"got <{type(point).__name__}>")

    if self._lbl.Lat.nbrows - 1 != len(point):
      raise ValueError(f"incompatible dimensions: {self._lbl.Lat.nbrows - 1} "
                       f"(LBL) vs {len(point)} (vector)")

    return self._lbl.contains_point(point)

  def zdomain(self):
    """
    Calculate the entire domain (Z-domain) of this LBL,
    by eliminating the existential variables.

    Returns:
      LBL: The corresponding Z-domain
    """
    return LBL(self._lbl.z_domain())

  def disjoint(self):
    """
    Calculate the disjoint union of this LBL, each single LBL of it is not
    intersecting the other ones.

    Returns:
      LBL: disjoint LBL union
    """
    return LBL(self._lbl.disjoint())

  def image(self, transfo):
    """
    Calculates the image of this LBL under an affine transformation.

    Args:
      transfo (Transfo): The transformation to apply

    Returns:
      LBL: The image of the LBL under the transformation

    Exemple:
      >>> a = LBLRead("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
      >>> f = Transfo("(i,j -> 2i+1, i+j)")
      >>> print(a.image(f))
    """
    if not isinstance(transfo, Transfo):
      raise TypeError(
        f"LBL.image() expected a Transfo, got <{type(transfo).__name__}>"
      )
    return LBL(self._lbl.image(transfo._mat))

  def preimage(self, transfo):
    """
    Calculates the preimage of this LBL under an affine transformation.

    Args:
      transfo (Transfo): The transformation to be inverted

    Returns:
      LBL: The preimage of the LBL under the transformation
    """
    if not isinstance(transfo, Transfo):
      raise TypeError(
        f"LBL.preimage() expected a Transfo, got <{type(transfo).__name__}>"
      )
    return LBL(self._lbl.preimage(transfo._mat))

  def plot(self, *args, **kwargs):
    """
    Plots the LBL using the plotting library in lbl_plot.py (pyvista)

    Optional arguments:
    - show_points (boolean): draw points that are part of the LBL
      (default: True)
    - subplot (boolean): set this to True if you want to plot another LBL
        in the same window -it will be rendered later in that case
        (default: False)
    - show_grid (string): 'custom', 'full box', 'pyvista' or None
        can be extended in future versions (default: 'custom')
    - color (string): color used to plot this object (default: random for each
        sub-LBL of the union)
    - supplementary arguments are transmitted to the plotter function
      (pyvista's plotter.show(**kwargs))
      Example usage: you can set option "interactive_update=True" to
      continue running python and get another plot window. Rendering will be
      ensured by the latest call to LBL.plot(). All currently opened plot
      windows will be closed when one of them is closed.
    """
    lbl_plot(self, *args, **kwargs)

  def sLBL_list(self):
    """
    Build a python list of single LBLs.

    lbl[i] is the i-th element of the union of single LBLs.
    """
    result = []

    node = self._lbl
    while node is not None:
      poly = node.P
      while poly is not None:
        result.append(LBL(pl.LBLAlloc(node.Lat, poly.copy_single_pol())))
        poly = poly.next
      node = node.next

    return result


  ######################## LBL methods overriding ########################
  def __iter__(self):
    """
    Lists all the integer points of the LBL.
    Iterates through the integer points of the polyhedron and calculates their
    image under the lattice.
    Throws a ValueError if the polyhedron is unbounded.

    Yields:
      tuple: Coordinates of the integer point in the lattice image
    """

    # iterate over single LBLs and call an iterator over each of them
    seen = set()
    node = self._lbl

    while node is not None:
      poly = node.P
      while poly is not None:
        yield from _iter_single_lbl(node.Lat, poly, seen)
        poly = poly.next
      node = node.next
  
  def __add__(self, other):
    """Union : a + b  ≡  a.union(b)"""
    return self.union(other)

  def __sub__(self, other):
    """Difference : a - b  ≡  a.difference(b)"""
    return self.difference(other)

  def __mul__(self, other):
    """Intersection : a * b  ≡  a.intersection(b)"""
    return self.intersection(other)

  def __eq__(self, other):
    """Equality : A == B  ≡  (A ⊆ B and B ⊆ A)"""
    if not isinstance(other, LBL):
      return NotImplemented
    return self.included(other) and other.included(self)

  def __contains__(self, other):
    """
    Inclusion test, geometric OR pythonic:
    - lbl1 in lbl2 : tests whether lbl1 ⊆ lbl2
    - point in lbl : tests whether a point (iterable) is in lbl
    """
    if isinstance(other, LBL):
      return other.included(self)
    else:
      return self.contains_point(other)

  def __bool__(self):
    # True if the LBL is not empty
    return self._lbl is not None and self._lbl.P != None


# ──────────────────────────────────────────────────────────────
# Transformation class
# ──────────────────────────────────────────────────────────────
class Transfo:
  """
  Linear affine transformation on integer points.

  A Transfo is defined by a symbolic notation of the form
  “(i,j -> 3i+1, 2i+5j)”, which is stored as a PolyLib matrix.

  Usage:
    f = Transfo("(i,j -> 3i+1, 2i+5j)")
    b = a.image(f)      # or f(a)         : image    of LBL a under f
    c = a.preimage(f)   # or f.preimage(a): preimage of LBL a under f
  """
  _mat:pl.Matrix = None
  def __init__(self, s=None):
    """
    Initializes a Transfo from a symbolic string.

    Args:
      s (str, optional): String in the format “(var1, var2, ... ->
                               expr1, expr2, ...)”.
                 If None, creates an empty Transfo.
    """
    if s is not None:
      self._mat = io.TransfoRead(s)

  def __repr__(self):
    """
    Returns the string representation of the transformation.

    Returns:
      str: example: “(i, j -> 3i+1, 2i+5j)”
    """
    if self._mat is None:
      return "None"
    return io.TransfoRepr(self._mat)

  def __call__(self, a):
    if not isinstance(a, LBL):
      raise NotImplemented
    return a.image(self)

  def preimage(self, a:LBL):
    """
    Calculates the preimage of LBL a under this affine transformation.

    Args:
      a (LBL): The source LBL

    Returns:
      LBL: The preimage of a under the transformation
    """
    if not isinstance(a, LBL):
      raise TypeError(
        f"Transfo.preimage() expected an LBL, got <{type(a).__name__}>"
      )
    return a.preimage(self)


  def compose(self, other):
    """
    Composes two transformations: self * other.

    Args:
      other (Transfo): The other transformation

    Returns:
      Transfo: The composition of the two transformations
    """
    if not isinstance(other, Transfo):
      raise TypeError(
        f"Transfo.compose() expected a Transfo, got <{type(other).__name__}>"
      )
    result = Transfo()
    result._mat = self._mat.multiply(other._mat)
    return result

  def inverse(self):
    """
    Calculates the inverse of this transformation.

    WARNING: this is an *integer* inverse transformation, it does not mean that
    t.inverse() . t == Id, but it is a constant times the identity matrix!
    Example:
      t = Transfo("(i, j -> 3i+1, i+j)")
      print(t.inverse())      # prints out: (i, j -> i-1, -i+3j+1)
      print(t * t.inverse())  # prints out: (i, j -> 3i, 3j)
    The image of an LBL under t is equal to its preimage under t.inverse();
    but the opposite is NOT true! If pA = (preimage of A under t) and iA =
    (image of A under t.inverse()), then the points of pA are a subset of the
    points of iA divided by the homothety ratio (t * t.inverse()).

    Returns:
      Transfo: The inverse transformation

    Raises:
      RuntimeError: If the matrix is not invertible
    """
    result = Transfo()
    result._mat = self._mat.inverse()
    return result

  def __mul__(self, other):
    """Composition : f * g  ≡  f.compose(g)"""
    return self.compose(other)


def _iter_single_lbl(lat, poly, seen):
  """
  Enumerates all the points of a single LBL (lat, poly)

  Yields:
    tuple: Coordinates of the point in the LBL
      (image by lat of an integer point of poly)

  Keeps the points that have been already seen once in the set 'seen'.
  """

  def get_bounds(scan_poly, k, point):
    cmat = scan_poly.constraint
    lo = None
    hi = None
    if not cmat:
      return lo, hi
    for r in range(scan_poly.nbconstraints):
      coeff_k = cmat[r, k + 1]
      if coeff_k == 0:
        continue
      val = cmat[r, dim + 1]
      for d in range(k):
        val += cmat[r, d + 1] * point[d]
      if coeff_k > 0:
        lo_c = math.ceil(-val / coeff_k)
        lo = lo_c if lo is None else max(lo, lo_c)
      else:
        hi_c = math.floor(val / (-coeff_k))
        hi = hi_c if hi is None else min(hi, hi_c)
    return lo, hi

  def _enumerate(k, point):
    if k == dim:
      # got a point being scanned
      img = []
      nb_vars = lat.nbcolumns - 1
      for r in range(lat.nbrows - 1):
        v = sum(lat[r, c] * point[c] for c in range(nb_vars))
        v += lat[r, nb_vars]
        img.append(v)
      pt = tuple(img)
      if pt not in seen:
        seen.add(pt)
        yield pt
      return
    if not scan_list:
      # got an early exit during plot (window closed?)
      return
    # continue scanning inner dimensions
    lo, hi = get_bounds(scan_list[k], k, point)
    if lo is None or hi is None:
      return
    for val in range(lo, hi + 1):
      yield from _enumerate(k + 1, point + [val])


  # _iter_single_lbl(self, lat, poly, seen) starts here:
  dim = poly.dimension

  if not poly.is_bounded():
    raise ValueError("unbounded; enumeration is impossible")

  # Use Polyhedron_Scan to get a scanning loop
  ctx_mat = pl.matrix_read_from_string("0 2\n")
  ctx = pl.constraints2polyhedron(ctx_mat)
  scan = poly.scan(ctx)

  # Collect scan polyhedra into a list (one per dimension)
  scan_list = []
  while scan is not None:
    scan_list.append(scan)
    scan = scan.next
  yield from _enumerate(0, [])
