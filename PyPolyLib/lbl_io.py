""" Various input/output functions for PyPolyLib.

The four functions (_LBL | _Transfo)(Read | Repr) are defined
"""


import re
import pypolylib_core as pl

def LBLRead(s):
    """
    Parse a symbolic string and constructs the internal LBL.

    Args:
        s (str): String in the format “{(expr1, ...) | constraints}”
                    Variables are lowercase letters (i, j, k, ...).
                    Supported constraints: <=, >=, =, <, >, double inequalities.

    Returns:
        pypolylib_core.LBL: The corresponding C LBL object
    """
    for i, c in enumerate(s):
        if ord(c) > 127:
            raise ValueError(f"non ASCII char {c!r} at position {i} "
                    f"(U+{ord(c):04X})\n")

    s = s.strip()
    if not (s.startswith('{') and s.endswith('}')):
        raise ValueError("The string must begin with ‘{’ and end with '}'")
    s = s[1:-1]

    if '|' not in s:
        raise ValueError("The separator ‘|’ is missing")
    lhs_str, rhs_str = s.split('|', 1)
    lhs_str = lhs_str.strip()
    rhs_str = rhs_str.strip()

    # enable multicharacter variables
    variables = re.findall(r'[a-zA-Z_][a-zA-Z0-9_]*(?![a-zA-Z0-9_])', s)
    # remove duplicates (but keep order)
    variables = list(dict.fromkeys(variables))
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

    # right hand side
    if rhs_str == "<empty>":
        return pl.LBLAlloc(lat, None)

    # Constraints Matrix
    constraint_rows = _parse_constraints(rhs_str, variables)
    nb_rows_poly = len(constraint_rows)
    nb_cols_poly = n + 2  # type + variables + constante

    poly_values = []
    for row in constraint_rows:
        poly_values.extend(row)

    cmat_str = _matrix_to_polylib_string(nb_rows_poly, nb_cols_poly, poly_values)
    cmat = pl.MatrixReadFromString(cmat_str)

    poly = pl.Constraints2Polyhedron(cmat)
    return pl.LBLAlloc(lat, poly)


def LBLRepr(lbl):
    """
    Constructs a string representation of an LBL.
    Ex: {(3i) | 0 <= i <= 5}

    Args:
        node: LBL node (pypolylib_core.LBL)
    """
    parts = []
    node = lbl
    while node is not None:
        P = node.P
        if P is None:
            return(_single_lbl_repr(node.Lat, P))
        while P:
            parts.append(_single_lbl_repr(node.Lat, P))
            P = P.next
        node = node.next
    return " UNION\n  " .join(parts)


def _single_lbl_repr(lat, poly):
    """
    Constructs the string representation of a single LBL node.

    Args:
        lat: lattice function (pypolylib_core.Matrix)
        poly: polyhedron (pypolylib_core.Polyhedron)

    Returns:
        str: Symbolic representation, e.g., “{(3i) | i >= 0, i <= 5}”
    """
    # lat = node.Lat
    # poly = node.P

    nb_out = lat.nbrows - 1      # number of outputs (excluding the homogeneous line)
    nb_vars = lat.nbcolumns - 1  # number of input variables

    # Variable Names : i, j, k, l, ...
    var_names = [chr(ord('i') + v) for v in range(nb_vars)]

    # ── Left side: output expressions ──
    out_exprs = []
    for r in range(nb_out):
        terms = []
        for c in range(nb_vars):
            coef = lat[r, c]
            if coef == 0:
                continue
            vname = var_names[c]
            if coef == 1:
                terms.append(vname)
            elif coef == -1:
                terms.append(f"-{vname}")
            else:
                terms.append(f"{coef}{vname}")
        cte = lat[r, nb_vars]
        if cte != 0:
            terms.append(str(cte))
        if not terms:
            out_exprs.append("0")
        else:
            expr = terms[0]
            for t in terms[1:]:
                if t.startswith('-'):
                    expr += t
                else:
                    expr += "+" + t
            out_exprs.append(expr)

    lhs = ", ".join(out_exprs)
    # if len(out_exprs) > 1 or True:
    lhs = f"({lhs})"

    # ── Right side: constraints ──
    if poly is None:
        # special string for empty LBLs: "{(_, _) | <empty>}" (2D example)
        return "{(" + ", ".join("_" for _ in range(nb_out)) + ") | <empty>}"

    nb_constraints = poly.nbconstraints
    dim = poly.dimension
    cmat = poly.constraint

    # Names of the polyhedron's variables
    poly_vars = [chr(ord('i') + v) for v in range(dim)]

    constraints = []
    for r in range(nb_constraints):
        eq_type = cmat[r, 0]  # 0 = equality, 1 = inequality
        coeffs = [cmat[r, c+1] for c in range(dim)]
        cte = cmat[r, dim+1]

        if coeffs == [0]*dim and cte == 1:
            # positivity constraint, don't print.
            continue

        # Construire l'expression : sum(coeff*var) + cte >= 0
        terms = []
        for c, coef in enumerate(coeffs):
            if coef == 0:
                continue
            vname = poly_vars[c]
            if coef == 1:
                terms.append(vname)
            elif coef == -1:
                terms.append(f"-{vname}")
            else:
                terms.append(f"{coef}{vname}")

        if eq_type == 0:
            # Equality : expr + cte = 0
            expr = _terms_to_str(terms)
            if cte != 0:
                constraints.append(f"{expr} = {-cte}")
            else:
                constraints.append(f"{expr} = 0")
        else:
            # Inequality: expr + cte >= 0
            # We want to write this in the form a <= var <= b if possible
            non_zero = [(c, coeffs[c]) for c in range(dim) if coeffs[c] != 0]
            if len(non_zero) == 1:
                c, coef = non_zero[0]
                vname = poly_vars[c]
                if coef > 0:
                    # coef*var + cte >= 0 -> var >= -cte/coef
                    bound = -cte // coef
                    constraints.append(f"{vname} >= {bound}")
                else:
                    # coef*var + cte >= 0 -> var <= cte/(-coef)
                    bound = cte // (-coef)
                    constraints.append(f"{vname} <= {bound}")
            else:
                expr = _terms_to_str(terms)
                if cte > 0:
                    constraints.append(f"{expr} + {cte} >= 0")
                elif cte < 0:
                    constraints.append(f"{expr} - {-cte} >= 0")
                else:
                    constraints.append(f"{expr} >= 0")

    rhs = ", ".join(constraints)
    return "{" + lhs + " | " + rhs + "}"


def _terms_to_str(terms):
    """
    Converts a list of terms into a linear expression.

    Args:
        terms (list): List of strings, e.g., [“2i”, “-3j”, “1”]

    Returns:
        str: Linear expression, e.g., “2i-3j+1”
    """
    if not terms:
        return "0"
    expr = terms[0]
    for t in terms[1:]:
        if t.startswith('-'):
            expr += t
        else:
            expr += "+" + t
    return expr


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
        str: String in PolyLib format, e.g., “2 3\n1 0 0\n0 1 0\n”
    """
    s = f"{nb_rows} {nb_cols}\n"
    idx = 0
    for i in range(nb_rows):
        s += " ".join(str(values[idx + j]) for j in range(nb_cols)) + "\n"
        idx += nb_cols
    return s

def TransfoRead(s):
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

def TransfoRepr(mat):

    nb_out = mat.nbrows - 1
    nb_vars = mat.nbcolumns - 1
    var_names = [chr(ord('i') + v) for v in range(nb_vars)]

    lhs = ", ".join(var_names)

    out_exprs = []
    for r in range(nb_out):
        terms = []
        for c in range(nb_vars):
            coef = mat[r, c]
            if coef == 0:
                continue
            vname = var_names[c]
            if coef == 1:    terms.append(vname)
            elif coef == -1: terms.append(f"-{vname}")
            else:            terms.append(f"{coef}{vname}")
        cte = mat[r, nb_vars]
        if cte != 0:
            terms.append(str(cte))
        out_exprs.append(_terms_to_str(terms) if terms else "0")

    rhs = ", ".join(out_exprs)
    return f"({lhs} -> {rhs})"
