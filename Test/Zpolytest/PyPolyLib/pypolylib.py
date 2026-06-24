"""
pypolylib.py — Python interface to the PolyLib library (Z-polyhedra).

This module exposes the LBL (Lattice-Based Lattice) and Transfo types
for manipulating sets of integer points (Z-polyhedra) in Python.

Typical usage:
    from pypolylib import LBLRead, Transfo

    #Create an LBL
    a = LBLRead("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
    b = LBLRead("{(2i, j) | 1 <= i <= 50, j = 10}")

    #  Set operations
    inter = a * b          # intersection
    diff  = a - b          # différence
    union = a + b          # union

    #  Check for inclusion
    print(a.included(b))   # True or False

    #Image under a transformation
    f = Transfo("(i,j -> 2i+1, i+3j)")
    print(a.image(f))      # image of a under f
    print(a.preimage(f))   # preimage of a under f

    # Composing and inverting transformations
    g = Transfo("(i,j -> i+1, 2j)")
    h = f * g              # composition
    fi = f.inverse()       # inverse

    # Iterate over integer points
    for pt in a:
        print(pt)          # prints each integer point

    # Plot (2D only)
    a.plot()
"""
#export LD_LIBRARY_PATH=$(pwd)/../../../.libs


import pypolylib_core as pl
import re
import math
from lbl_repr import _lbl_repr
from lbl_repr import _terms_to_str
from lbl_plot import lbl_plot


# ──────────────────────────────────────────────────────────────
# Internal parsing functions
# ──────────────────────────────────────────────────────────────

def _parse_linear(expr, variables):
    """
    Parse a linear expression in terms of the given variables.

    Args:
        expr (str): Linear expression, e.g.: “2i+3j-1”
        variables (list): List of variable names, e.g.: [‘i’, ‘j’]

    Returns:
        dict: Coefficients per variable + constant, e.g.: {‘i’:2, ‘j’:3, ‘cte’:-1}
    """
    coeffs = {v: 0 for v in variables}
    coeffs['cte'] = 0
    expr = expr.replace(' ', '').replace('-', '+-')
    terms = [t for t in expr.split('+') if t != '']
    for term in terms:
        if term == '' or term == '-':
            continue
        matched = False
        for v in variables:
            m = re.fullmatch(r'([+-]?\d*)' + re.escape(v), term)
            if m:
                c = m.group(1)
                if c in ('', '+'):   c = 1
                elif c == '-':       c = -1
                else:                c = int(c)
                coeffs[v] += c
                matched = True
                break
        if not matched:
            try:
                coeffs['cte'] += int(term)
            except ValueError:
                raise ValueError(f"Misunderstood Term : '{term}' in '{expr}'")
    return coeffs


def _parse_lhs(lhs_str, variables):
    """
    Parses the left-hand side of an LBL: “(expr1, expr2, ...)”.

    Args:
        lhs_str (str): Left-hand side, e.g.: “(i+j)” or “(2i, j+1)”
        variables (list): List of variable names

    Returns:
        list: List of dictionaries containing coefficients, one per output expression
    """
    lhs_str = lhs_str.strip()
    if lhs_str.startswith('(') and lhs_str.endswith(')'):
        lhs_str = lhs_str[1:-1]
    exprs = [e.strip() for e in lhs_str.split(',')]
    return [_parse_linear(e, variables) for e in exprs]


def _parse_constraints(rhs_str, variables):
    """
    Parse the constraints of an LBL into lines in PolyLib format.

    Supports: >=, <=, =, >, <, double inequalities (e.g., 0 <= i <= 10).
    Each returned line has the form [type, coeff_v1, ..., coeff_vn, constant]
    where type=0 for equality and type=1 for inequality (>= 0).

    Args:
        rhs_str (str): Right-hand side, e.g., “0 <= i <= 10, i = 2j”
        variables (list): List of variable names

    Returns:
        list: List of lists of integers representing the PolyLib constraints
    """
    rows = []
    constraints = [c.strip() for c in rhs_str.split(',')]
    for c in constraints:
        if not c:
            continue
        if re.search(r'(?<![<>])=(?!=)', c) and '<=' not in c and '>=' not in c:
            parts = re.split(r'=', c)
            lhs = _parse_linear(parts[0].strip(), variables)
            rhs = _parse_linear(parts[1].strip(), variables)
            row = [0]
            for v in variables:
                row.append(lhs[v] - rhs[v])
            row.append(lhs['cte'] - rhs['cte'])
            rows.append(row)
        elif c.count('<=') == 2:
            parts = re.split(r'<=', c)
            lhs0 = _parse_linear(parts[0].strip(), variables)
            mid  = _parse_linear(parts[1].strip(), variables)
            rhs0 = _parse_linear(parts[2].strip(), variables)
            row1 = [1] + [mid[v] - lhs0[v] for v in variables] + [mid['cte'] - lhs0['cte']]
            row2 = [1] + [rhs0[v] - mid[v] for v in variables] + [rhs0['cte'] - mid['cte']]
            rows.append(row1)
            rows.append(row2)
        elif c.count('>=') == 2:
            parts = re.split(r'>=', c)
            lhs0 = _parse_linear(parts[0].strip(), variables)
            mid  = _parse_linear(parts[1].strip(), variables)
            rhs0 = _parse_linear(parts[2].strip(), variables)
            row1 = [1] + [lhs0[v] - mid[v] for v in variables] + [lhs0['cte'] - mid['cte']]
            row2 = [1] + [mid[v] - rhs0[v] for v in variables] + [mid['cte'] - rhs0['cte']]
            rows.append(row1)
            rows.append(row2)
        elif '<=' in c:
            parts = re.split(r'<=', c)
            lhs = _parse_linear(parts[0].strip(), variables)
            rhs = _parse_linear(parts[1].strip(), variables)
            row = [1] + [rhs[v] - lhs[v] for v in variables] + [rhs['cte'] - lhs['cte']]
            rows.append(row)
        elif '>=' in c:
            parts = re.split(r'>=', c)
            lhs = _parse_linear(parts[0].strip(), variables)
            rhs = _parse_linear(parts[1].strip(), variables)
            row = [1] + [lhs[v] - rhs[v] for v in variables] + [lhs['cte'] - rhs['cte']]
            rows.append(row)
        elif c.count('<') == 2 and '<=' not in c:
            parts = re.split(r'<', c)
            lhs0 = _parse_linear(parts[0].strip(), variables)
            mid  = _parse_linear(parts[1].strip(), variables)
            rhs0 = _parse_linear(parts[2].strip(), variables)
            row1 = [1] + [mid[v] - lhs0[v] for v in variables] + [mid['cte'] - lhs0['cte'] - 1]
            row2 = [1] + [rhs0[v] - mid[v] for v in variables] + [rhs0['cte'] - mid['cte'] - 1]
            rows.append(row1)
            rows.append(row2)
        elif c.count('>') == 2 and '>=' not in c:
            parts = re.split(r'>', c)
            lhs0 = _parse_linear(parts[0].strip(), variables)
            mid  = _parse_linear(parts[1].strip(), variables)
            rhs0 = _parse_linear(parts[2].strip(), variables)
            row1 = [1] + [lhs0[v] - mid[v] for v in variables] + [lhs0['cte'] - mid['cte'] - 1]
            row2 = [1] + [mid[v] - rhs0[v] for v in variables] + [mid['cte'] - rhs0['cte'] - 1]
            rows.append(row1)
            rows.append(row2)
        elif '<' in c and '<=' not in c:
            parts = re.split(r'<', c)
            lhs = _parse_linear(parts[0].strip(), variables)
            rhs = _parse_linear(parts[1].strip(), variables)
            row = [1] + [rhs[v] - lhs[v] for v in variables] + [rhs['cte'] - lhs['cte'] - 1]
            rows.append(row)
        elif '>' in c and '>=' not in c:
            parts = re.split(r'>', c)
            lhs = _parse_linear(parts[0].strip(), variables)
            rhs = _parse_linear(parts[1].strip(), variables)
            row = [1] + [lhs[v] - rhs[v] for v in variables] + [lhs['cte'] - rhs['cte'] - 1]
            rows.append(row)
        else:
            raise ValueError(f"Misunderstood Constraint: '{c}'")
    return rows


def _matrix_to_polylib_string(nb_rows, nb_cols, values):
    """
    Constructs a string in PolyLib format for Matrix_Read.

    Args:
        nb_rows (int): Number of rows
        nb_cols (int): Number of columns
        values (list): Matrix values, row by row

    Returns:
        str: String in PolyLib format, e.g., “2 3\\n1 0 0\\n0 1 0\\n”
    """
    s = f"{nb_rows} {nb_cols}\n"
    idx = 0
    for i in range(nb_rows):
        s += " ".join(str(values[idx + j]) for j in range(nb_cols)) + "\n"
        idx += nb_cols
    return s


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
    Z-polyhedron (Lattice-Based Lattice) — set of integer points.

    An LBL is defined by a lattice (transformation matrix) and
    a constraint polyhedron. It represents the image of a set
    of integer points under a linear affine function.

    Creation :
        a = LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
        # or :
        a = LBLRead("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")

    Operations :
        a + b   # union
        a - b   # difference
        a * b   # intersection
    """

    def _LBLRead(self, s):
        """
        Parse a symbolic string and constructs the internal LBL.

        Args:
            s (str): String in the format “{(expr1, ...) | constraints}”
                     Variables are lowercase letters (i, j, k, ...).
                     Supported constraints: <=, >=, =, <, >, double inequalities.

        Returns:
            pypolylib_core.LBL: The corresponding C LBL object
        """
        s = s.strip()
        if not (s.startswith('{') and s.endswith('}')):
            raise ValueError("The string must begin with ‘{’ and end with '}'")
        s = s[1:-1]

        if '|' not in s:
            raise ValueError("The separator ‘|’ is missing")
        lhs_str, rhs_str = s.split('|', 1)
        lhs_str = lhs_str.strip()
        rhs_str = rhs_str.strip()

        all_text = lhs_str + rhs_str
        variables = sorted(set(re.findall(r'[a-z](?![a-z])', all_text)))
        n = len(variables)

        # Lattice matrix
        lhs_coeffs = _parse_lhs(lhs_str, variables)
        nb_out = len(lhs_coeffs)
        nb_rows_lat = nb_out + 1
        nb_cols_lat = n + 1

        lat_values = []
        for coeffs in lhs_coeffs:
            for v in variables:
                lat_values.append(coeffs[v])
            lat_values.append(coeffs['cte'])
        # uniform line
        for j in range(n):
            lat_values.append(0)
        lat_values.append(1)

        lat_str = _matrix_to_polylib_string(nb_rows_lat, nb_cols_lat, lat_values)
        lat = pl.MatrixReadFromString(lat_str)

        # Constraints Matrix
        constraint_rows = _parse_constraints(rhs_str, variables)
        nb_rows_poly = len(constraint_rows)
        nb_cols_poly = n + 2  # type + variables + constante

        poly_values = []
        for row in constraint_rows:
            poly_values.extend(row)

        cmat_str = _matrix_to_polylib_string(nb_rows_poly, nb_cols_poly, poly_values)
        cmat = pl.MatrixReadFromString(cmat_str)

        poly = pl.Constraints2Polyhedron(cmat, 0)
        return pl.LBLAlloc(lat, poly)

    def __init__(self, s=None):
        """
        Initializes an LBL, optionally from a symbolic string.

        Arguments:
            s (str, optional): string in the format “{(expr) | constraints}”.
                               If None, creates an empty LBL.
        """
        self._lbl = None
        if s is not None:
            self._lbl = self._LBLRead(s)

    def Print(self):
        """Displays the LBL in symbolic form."""
        print(self.__repr__())

    def __repr__(self):
        """Returns a symbolic representation of the LBL, e.g., {(3i) | i >= 0, i <= 5}"""
        if self._lbl is None:
            return "LBL(empty)"
        return _lbl_repr(self._lbl)

    def intersection(self, other):
        """
        Calculates the intersection of this LBL with another one.

        Args:
            other (LBL): The other LBL (same dimension)

        Returns:
            LBL: The intersection of the two LBLs
        """
        if self._lbl.P.dimension != other._lbl.P.dimension:
            raise ValueError(f"Incompatible dimensions: {self._lbl.P.dimension} vs {other._lbl.P.dimension}")
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
        if self._lbl.P.dimension != other._lbl.P.dimension:
            raise ValueError(f"Incompatible dimensions: {self._lbl.P.dimension} vs {other._lbl.P.dimension}")
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
        if self._lbl.P.dimension != other._lbl.P.dimension:
            raise ValueError(f"Incompatible dimensions: {self._lbl.P.dimension} vs {other._lbl.P.dimension}")
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

    def plot(self):
        """Displays the LBL using matplotlib (2D only)."""
        lbl_plot(self)

    def __iter__(self):
        """
        Lists all the integer points of the LBL.
        Iterates through the integer points of the polyhedron and calculates their image under the lattice.
        Throws an error if the polyhedron is unbounded.

        Yields:
            tuple: Coordinates of the integer point in the lattice image
        """
        seen = set()
        node = self._lbl

        while node is not None:
            lat = node.Lat
            poly = node.P

            # Boucle interne : union de polyèdres dans le même node
            while poly is not None:
                nb_out = lat.nbrows - 1
                nb_vars = lat.nbcolumns - 1
                dim = poly.dimension

                # Check whether the polyhedron is bounded
                rmat = poly.rays
                for r in range(poly.nbrays):
                    last_col = pl.MatrixGetValue(rmat, r, poly.dimension + 1)
                    if last_col == 0:  # rayon infini
                        raise ValueError("The polyhedron is unbounded; enumeration is impossible")

                # # Extract the constraints: coefficients * x + constant >= 0
                # constraints = []
                # nb_constraints = poly.nbconstraints
                # cmat = poly.constraints
                # for r in range(nb_constraints):
                #     eq_type = pl.MatrixGetValue(cmat, r, 0)
                #     coeffs = [pl.MatrixGetValue(cmat, r, c+1) for c in range(dim)]
                #     cte = pl.MatrixGetValue(cmat, r, dim+1)
                #     constraints.append((eq_type, coeffs, cte))

                    # Use Polyhedron_Scan for optimal per-dimension bounds
                ctx_mat = pl.MatrixReadFromString("0 2\n")
                ctx = pl.Constraints2Polyhedron(ctx_mat)
                poly_single = pl.PolyhedronCopy(poly)
                scan = pl.PolyhedronScan(poly_single, ctx, 1024)

                # Collect scan polyhedra into a list (one per dimension)
                scan_list = []
                s = scan
                while s is not None:
                    scan_list.append(s)
                    s = s.next

                def get_bounds(scan_poly, k, point):
                    cmat = scan_poly.constraints
                    lo = None
                    hi = None
                    for r in range(scan_poly.nbconstraints):
                        coeff_k = pl.MatrixGetValue(cmat, r, k + 1)
                        if coeff_k == 0:
                            continue
                        val = pl.MatrixGetValue(cmat, r, dim + 1)
                        for d in range(k):
                            val += pl.MatrixGetValue(cmat, r, d + 1) * point[d]
                        if coeff_k > 0:
                            lo_c = math.ceil(-val / coeff_k)
                            lo = lo_c if lo is None else max(lo, lo_c)
                        else:
                            hi_c = math.floor(val / (-coeff_k))
                            hi = hi_c if hi is None else min(hi, hi_c)
                    return lo, hi

                def _enumerate(k, point):
                    if k == dim:
                        img = []
                        for r in range(nb_out):
                            v = sum(pl.MatrixGetValue(lat, r, c) * point[c] for c in range(nb_vars))
                            v += pl.MatrixGetValue(lat, r, nb_vars)
                            img.append(v)
                        pt = tuple(img)
                        if pt not in seen:
                            seen.add(pt)
                            yield pt
                        return
                    lo, hi = get_bounds(scan_list[k], k, point)
                    if lo is None or hi is None:
                        return
                    for val in range(lo, hi + 1):
                        yield from _enumerate(k + 1, point + [val])

                yield from _enumerate(0, [])

                poly = poly.next  # polyèdre suivant dans le même node
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
        f = Transfo(“(i,j -> 3i+1, 2i+5j)”)
        b = a.image(f)      # image of an LBL under f
        c = a.preimage(f)   # preimage of an LBL under f
    """

    def __init__(self, s=None):
        """
        Initializes a Transfo from a symbolic string.

        Args:
            s (str, optional): String in the format “(var1, var2, ... -> expr1, expr2, ...)”.
                               If None, creates an empty Transfo.
        """
        self._mat = None
        if s is not None:
            self._mat = self._parse(s)

    def _parse(self, s):
        """
        Parse the symbolic string and constructs the internal matrix.

        Args:
            s (str): String in the format “(i,j -> 3i+1, 2i+5j)”

        Returns:
            pypolylib_core.Matrix: The PolyLib transformation matrix
        """
        s = s.strip()
        if s.startswith('(') and s.endswith(')'):
            s = s[1:-1]
        if '->' not in s:
            raise ValueError("'->' is missing from the transformation")

        lhs_str, rhs_str = s.split('->', 1)

        # Input variables : i, j, k, ...
        variables = [v.strip() for v in lhs_str.split(',')]
        n = len(variables)

        # Output expressions
        out_exprs = [e.strip() for e in rhs_str.split(',')]
        nb_out = len(out_exprs)

        # Matrix (nb_out+1) x (n+1)
        nb_rows = nb_out + 1
        nb_cols = n + 1

        values = []
        for expr in out_exprs:
            coeffs = _parse_linear(expr, variables)
            for v in variables:
                values.append(coeffs[v])
            values.append(coeffs['cte'])
        # uniform lign
        for j in range(n):
            values.append(0)
        values.append(1)

        mat_str = _matrix_to_polylib_string(nb_rows, nb_cols, values)
        return pl.MatrixReadFromString(mat_str)

    def __repr__(self):
        """
        Returns the symbolic representation of the transformation.

        Returns:
            str: Example: “(i, j -> 3i+1, 2i+5j)”
        """
        if self._mat is None:
            return "Transfo(vide)"
        nb_out = self._mat.nbrows - 1
        nb_vars = self._mat.nbcolumns - 1
        var_names = [chr(ord('i') + v) for v in range(nb_vars)]

        lhs = ", ".join(var_names)

        out_exprs = []
        for r in range(nb_out):
            terms = []
            for c in range(nb_vars):
                coef = pl.MatrixGetValue(self._mat, r, c)
                if coef == 0:
                    continue
                vname = var_names[c]
                if coef == 1:    terms.append(vname)
                elif coef == -1: terms.append(f"-{vname}")
                else:            terms.append(f"{coef}{vname}")
            cte = pl.MatrixGetValue(self._mat, r, nb_vars)
            if cte != 0:
                terms.append(str(cte))
            out_exprs.append(_terms_to_str(terms) if terms else "0")

        rhs = ", ".join(out_exprs)
        return f"({lhs} -> {rhs})"

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