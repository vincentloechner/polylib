import pypolylib
import re


# Utilitaire : parser une expression lineaire

def _parse_linear(expr, variables):
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
                raise ValueError(f"Terme incompris : '{term}' dans '{expr}'")
    return coeffs


def _parse_lhs(lhs_str, variables):
    lhs_str = lhs_str.strip()
    if lhs_str.startswith('(') and lhs_str.endswith(')'):
        lhs_str = lhs_str[1:-1]
    exprs = [e.strip() for e in lhs_str.split(',')]
    return [_parse_linear(e, variables) for e in exprs]


def _parse_constraints(rhs_str, variables):
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
        else:
            raise ValueError(f"Contrainte incomprise : '{c}'")
    return rows


def _matrix_to_polylib_string(nb_rows, nb_cols, values):
    """Construit une chaîne au format PolyLib pour Matrix_Read."""
    s = f"{nb_rows} {nb_cols}\n"
    idx = 0
    for i in range(nb_rows):
        s += " ".join(str(values[idx + j]) for j in range(nb_cols)) + "\n"
        idx += nb_cols
    return s


def LBLRead(s):
    """
    Parse une chaîne du type "{(2i, j) | 0 <= i <= 10, i = 2j}"
    et retourne un objet LBL de pypolylib.
    """
    s = s.strip()
    if not (s.startswith('{') and s.endswith('}')):
        raise ValueError("La chaîne doit commencer par '{' et finir par '}'")
    s = s[1:-1]

    if '|' not in s:
        raise ValueError("Il manque le séparateur '|'")
    lhs_str, rhs_str = s.split('|', 1)
    lhs_str = lhs_str.strip()
    rhs_str = rhs_str.strip()

    all_text = lhs_str + rhs_str
    variables = sorted(set(re.findall(r'[a-z](?![a-z])', all_text)))
    n = len(variables)

    #  Matrice lattice 
    lhs_coeffs = _parse_lhs(lhs_str, variables)
    nb_out = len(lhs_coeffs)
    nb_rows_lat = nb_out + 1
    nb_cols_lat = n + 1

    lat_values = []
    for coeffs in lhs_coeffs:
        for v in variables:
            lat_values.append(coeffs[v])
        lat_values.append(coeffs['cte'])
    # ligne homogene
    for j in range(n):
        lat_values.append(0)
    lat_values.append(1)

    lat_str = _matrix_to_polylib_string(nb_rows_lat, nb_cols_lat, lat_values)
    lat = pypolylib.MatrixReadFromString(lat_str)

    #  Matrice de contraintes 
    constraint_rows = _parse_constraints(rhs_str, variables)
    nb_rows_poly = len(constraint_rows)
    nb_cols_poly = n + 2  # type + variables + constante

    poly_values = []
    for row in constraint_rows:
        poly_values.extend(row)

    cmat_str = _matrix_to_polylib_string(nb_rows_poly, nb_cols_poly, poly_values)
    cmat = pypolylib.MatrixReadFromString(cmat_str)

    poly = pypolylib.Constraints2Polyhedron(cmat, 0)
    return pypolylib.LBLAlloc(lat, poly)


class LBL:
    """Classe Python enveloppant un LBL de PolyLib."""

    def __init__(self, s=None):
        self._lbl = None
        if s is not None:
            self._lbl = LBLRead(s)

    def Print(self):
        if self._lbl is None:
            print("LBL vide")
            return
        pypolylib.LBLPrint(self._lbl)

    def __repr__(self):
        if self._lbl is None:
            return "LBL(vide)"
        return "LBL(non vide)"