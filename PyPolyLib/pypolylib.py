"""
pypolylib.py — Python interface to the PolyLib library (LBLs).

This module exposes the LBL (image of Z-polyhedron) and Transfo (matrix)
types for manipulating sets of integer points in Python.

Typical usage:

```
from pypolylib import LBL, Transfo

# Create an LBL
a = LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
b = LBL("{(2i, j) | 1 <= i <= 50, j = 10}")

#  Set operations
inter = a * b          # intersection
diff  = a - b          # difference
union = a + b          # union

#  Check for inclusion
print(a.included(b))   # True or False

#  Check for equality (mutual inclusion)
print(a == b)   # True or False

# Image under a transformation
f = Transfo("(i,j -> 2i+1, i+3j)")
print(a.image(f))      # image of a under f
print(a.preimage(f))   # preimage of a under f

# Composing and inverting transformations
g = Transfo("(i,j -> i+1, 2j)")
h = f * g              # composition
fi = f.inverse()       # inverse

# Iterate over LBL points (bounded LBL)
for pt in a:
    print(pt)          # print each point (as a tuple of integers)
# get (python) set of points (bounded LBL)
s = set(a)

# Plot (2D and 3D only, bounded LBL only)
a.plot()
```
"""


import pypolylib_core as pl
import math
# from lbl_io import _lbl_repr, _terms_to_str
import lbl_io as io
from lbl_plot import lbl_plot



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


# ──────────────────────────────────────────────────────────────
# Classe LBL
# ──────────────────────────────────────────────────────────────

class LBL:
    """
    An LBL is a union of affine image of Z-domains — set of integer points.

    An LBL is stored as a polylib _lbl: a list of (lattice, polyhedral domain)

    Creation :
        a = LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
        # or :
        a = LBLRead("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")

    Operations :
        a + b   # union
        a - b   # difference
        a * b   # intersection
    """
    def __init__(self, s=None):
        """
        Initializes an LBL, optionally from a symbolic string.

        Arguments:
            s (str, optional): string in the format “{(expr) | constraints}”.
                               If None, creates an empty LBL.
        """
        self._lbl = None
        if s is not None:
            self._lbl = io._LBLRead(s)

    def print(self):
        """Displays the LBL in symbolic form."""
        print(self.__repr__())

    def __repr__(self):
        """Returns a symbolic representation of the LBL, e.g., {(3i) | i >= 0, i <= 5}"""
        if self._lbl is None:
            return "None"
        return io._LBLRepr(self._lbl)

    def intersection(self, other):
        """
        Calculates the intersection of this LBL with another one.

        Args:
            other (LBL): The other LBL (same dimension)

        Returns:
            LBL: The intersection of the two LBLs
        """
        if self._lbl.Lat.nbrows - 1 != other._lbl.Lat.nbrows - 1:
            raise ValueError(f"Incompatible dimensions: {self._lbl.Lat.nbrows - 1} vs {other._lbl.Lat.nbrows - 1}")
        result = LBL()
        result._lbl = self._lbl.intersection(other._lbl)
        return result

    def difference(self, other):
        """
        Calculates the difference between this LBL and another one (self - other).

        Args:
            other (LBL): The other LBL (same dimension)

        Returns:
            LBL: The difference between the two LBLs
        """
        if self._lbl.Lat.nbrows - 1 != other._lbl.Lat.nbrows - 1:
            raise ValueError(f"Incompatible dimensions: {self._lbl.Lat.nbrows - 1} vs {other._lbl.Lat.nbrows - 1}")
        result = LBL()
        result._lbl = self._lbl.difference(other._lbl)
        return result

    def union(self, other):
        """
        Calculates the union of this LBL with another one.

        Args:
            other (LBL): The other LBL (same dimension)

        Returns:
            LBL: The union of the two LBLs
        """
        if self._lbl.Lat.nbrows - 1 != other._lbl.Lat.nbrows - 1:
            raise ValueError(f"Incompatible dimensions: {self._lbl.Lat.nbrows - 1} vs {other._lbl.Lat.nbrows - 1}")
        result = LBL()
        result._lbl = self._lbl.union(other._lbl)
        return result

    def included(self, other):
        """
        Checks whether this LBL is contained within another one.

        Args:
            other (LBL): The other LBL

        Returns:
            bool: True if self ⊆ other, False otherwise
        """
        return self._lbl.included(other._lbl)

    def zdomain(self):
        """
        Calculate the entire domain (Z-domain) of this LBL,
        by eliminating the existential variables.

        Returns:
            LBL: The corresponding Z-domain
        """
        result = LBL()
        result._lbl = self._lbl.zdomain()
        return result

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
        result = LBL()
        result._lbl = pl.LBLImage(self._lbl, transfo._mat)
        return result

    def preimage(self, transfo):
        """
        Calculates the preimage of this LBL under an affine transformation.

        Args:
            transfo (Transfo): The transformation to be inverted

        Returns:
            LBL: The preimage of the LBL under the transformation
        """
        result = LBL()
        result._lbl = pl.LBLPreimage(self._lbl, transfo._mat)
        return result

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
        """Equality : A == B  ≡  A ⊆ B et B ⊆ A"""
        if not isinstance(other, LBL):
            return NotImplemented
        return self.included(other) and other.included(self)

    def __contains__(self, other):
        """Inclusion : other in self  ≡  other ⊆ self"""
        if not isinstance(other, LBL):
            return NotImplemented
        return other.included(self)

    def plot(self, *args, **kwargs):
        """Displays the LBL.
        
            Optionnal arguments:
            - show_points (boolean): draw points that are part of the LBL
                (default True)
            - subplot (boolean): set this to True if you want to plot
                another LBL in the same window (it will be rendered later in
                that case)
            - supplementary arguments are transmitted to the pyvista
                plotter.show() function. Example: you can set option
                 "interactive_update=True" to continue running python and get
                  another plot window activated.
        """
        lbl_plot(self, *args, **kwargs)

    def _iter_single_lbl(self, lat, poly, seen):
        """
        Lists all the integer points of a single LBL (lat, poly)
        Yields:
            tuple: Coordinates of the integer point in the lattice image
        """

        def get_bounds(scan_poly, k, point):
            cmat = scan_poly.constraint
            lo = None
            hi = None
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
            # else, continue scanning inside dimensions
            lo, hi = get_bounds(scan_list[k], k, point)
            if lo is None or hi is None:
                return
            for val in range(lo, hi + 1):
                yield from _enumerate(k + 1, point + [val])

        dim = poly.dimension

        # Check whether the polyhedron is bounded
        if poly.nbbid != 0:
            raise ValueError("unbounded; enumeration is impossible")
        if poly.ray[0, dim + 1] == 0:
            # infinite ray
            raise ValueError("unbounded; enumeration is impossible")

        # Use Polyhedron_Scan for optimal per-dimension bounds
        ctx_mat = pl.MatrixReadFromString("0 2\n")
        ctx = pl.Constraints2Polyhedron(ctx_mat)
        scan = pl.PolyhedronScan(poly, ctx, 1024)

        # Collect scan polyhedra into a list (one per dimension)
        scan_list = []
        while scan is not None:
            scan_list.append(scan)
            scan = scan.next
        yield from _enumerate(0, [])

    def __iter__(self):
        """
        Lists all the integer points of the LBL.
        Iterates through the integer points of the polyhedron and calculates their image under the lattice.
        Throws an error if the polyhedron is unbounded.

        Yields:
            tuple: Coordinates of the integer point in the lattice image
        """

        # iterate over single LBLs and call an iterator over each of them
        seen = set()
        node = self._lbl

        while node is not None:
            poly = node.P
            while poly is not None:
                yield from self._iter_single_lbl(node.Lat, poly, seen)

                poly = poly.next

            node = node.next


# ──────────────────────────────────────────────────────────────
# Classe Transfo
# ──────────────────────────────────────────────────────────────

class Transfo:
    """
    Linear affine transformation on integer points.

    A Transfo is defined by a symbolic notation of the form
    “(i,j -> 3i+1, 2i+5j)”, which is converted to a PolyLib matrix.

    Usage:
        f = Transfo("(i,j -> 3i+1, 2i+5j)")
        b = a.image(f)      # image of an LBL under f
        c = a.preimage(f)   # preimage of an LBL under f
    """
    def __init__(self, s=None):
        """
        Initializes a Transfo from a symbolic string.

        Args:
            s (str, optional): String in the format “(var1, var2, ... ->
                                                           expr1, expr2, ...)”.
                               If None, creates an empty Transfo.
        """
        self._mat = None
        if s is not None:
            self._mat = io._TransfoRead(s)

    def __repr__(self):
        """
        Returns the symbolic representation of the transformation.

        Returns:
            str: Example: “(i, j -> 3i+1, 2i+5j)”
        """
        if self._mat is None:
            return "None"
        return io._TransfoRepr(self._mat)

    def __call__(self, a):
        if not isinstance(a, LBL):
            raise TypeError(
                f"'{type(a).__name__}' object cannot be interpreted as an LBL")
        return a.image(self)

    def compose(self, other):
        """
        Composes two transformations: self * other.

        Args:
            other (Transfo): The other transformation

        Returns:
            Transfo: The composition of the two transformations
        """
        result = Transfo()
        result._mat = pl.MatrixProduct(self._mat, other._mat)
        return result

    def inverse(self):
        """
        Calculates the inverse of this transformation.

        Returns:
            Transfo: The inverse transformation

        Raises:
            RuntimeError: If the matrix is not invertible
        """
        result = Transfo()
        result._mat = pl.MatrixInverse(self._mat)
        return result

    def __mul__(self, other):
        """Composition : f * g  ≡  f.compose(g)"""
        return self.compose(other)