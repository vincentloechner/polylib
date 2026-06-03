import pypolylib
import re

# ──────────────────────────────────────────────
# Utilitaire : parser une expression linéaire
# ex: "2i+3j-1" -> {i:2, j:3, cte:-1}
# ──────────────────────────────────────────────
def _parse_linear(expr, variables):
    """
    Parse une expression linéaire en termes des variables données.
    Retourne un dict {var: coeff, ..., 'cte': valeur_constante}.
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
            # Cherche coeff*v, ex: "2i", "-3j", "i", "-i"
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
            # terme constant
            try:
                coeffs['cte'] += int(term)
            except ValueError:
                raise ValueError(f"Terme incompris : '{term}' dans '{expr}'")

    return coeffs


# ──────────────────────────────────────────────
# Parser la partie gauche : "(i+j)" ou "(2i, 0, 52i+66j)"
# Retourne une liste d'expressions linéaires (une par sortie)
# ──────────────────────────────────────────────
def _parse_lhs(lhs_str, variables):
    """
    Parse la partie gauche d'un LBL : "(expr1, expr2, ...)"
    Retourne une liste de dicts de coefficients.
    """
    lhs_str = lhs_str.strip()
    # Enlever les parenthèses extérieures
    if lhs_str.startswith('(') and lhs_str.endswith(')'):
        lhs_str = lhs_str[1:-1]
    exprs = [e.strip() for e in lhs_str.split(',')]
    return [_parse_linear(e, variables) for e in exprs]


# ──────────────────────────────────────────────
# Parser la partie droite : contraintes
# ex: "0 <= i <= 10, i = 2j"
# Retourne une liste de (type, coeffs) où type=0 égalité, 1 inégalité
# Format PolyLib : [type, coeff_v1, coeff_v2, ..., cte]
#   inégalité : expr >= 0
#   égalité   : expr = 0
# ──────────────────────────────────────────────
def _parse_constraints(rhs_str, variables):
    """
    Parse les contraintes et retourne une liste de lignes PolyLib.
    Chaque ligne = [eq_type, coeff_var1, ..., coeff_varn, constante]
    eq_type : 1 = inégalité (>=0), 0 = égalité (=0)
    """
    rows = []
    constraints = [c.strip() for c in rhs_str.split(',')]

    for c in constraints:
        if not c:
            continue

        # Égalité : "i = 2j" ou "i = 0"
        if re.search(r'(?<![<>])=(?!=)', c) and '<=' not in c and '>=' not in c:
            parts = re.split(r'=', c)
            lhs = _parse_linear(parts[0].strip(), variables)
            rhs = _parse_linear(parts[1].strip(), variables)
            # lhs - rhs = 0
            row = [0]
            for v in variables:
                row.append(lhs[v] - rhs[v])
            row.append(lhs['cte'] - rhs['cte'])
            rows.append(row)

        # Double inégalité : "0 <= i <= 10"
        elif c.count('<=') == 2 or c.count('>=') == 2:
            if '<=' in c:
                parts = re.split(r'<=', c)
                # parts[0] <= parts[1] <= parts[2]
                # => parts[1] - parts[0] >= 0
                lhs0 = _parse_linear(parts[0].strip(), variables)
                mid  = _parse_linear(parts[1].strip(), variables)
                rhs0 = _parse_linear(parts[2].strip(), variables)
                # mid - lhs0 >= 0
                row1 = [1]
                for v in variables:
                    row1.append(mid[v] - lhs0[v])
                row1.append(mid['cte'] - lhs0['cte'])
                rows.append(row1)
                # rhs0 - mid >= 0
                row2 = [1]
                for v in variables:
                    row2.append(rhs0[v] - mid[v])
                row2.append(rhs0['cte'] - mid['cte'])
                rows.append(row2)
            else:
                parts = re.split(r'>=', c)
                lhs0 = _parse_linear(parts[0].strip(), variables)
                mid  = _parse_linear(parts[1].strip(), variables)
                rhs0 = _parse_linear(parts[2].strip(), variables)
                row1 = [1]
                for v in variables:
                    row1.append(lhs0[v] - mid[v])
                row1.append(lhs0['cte'] - mid['cte'])
                rows.append(row1)
                row2 = [1]
                for v in variables:
                    row2.append(mid[v] - rhs0[v])
                row2.append(mid['cte'] - rhs0['cte'])
                rows.append(row2)

        # Inégalité simple : "i >= 0", "i <= 10"
        elif '<=' in c:
            parts = re.split(r'<=', c)
            lhs = _parse_linear(parts[0].strip(), variables)
            rhs = _parse_linear(parts[1].strip(), variables)
            # rhs - lhs >= 0
            row = [1]
            for v in variables:
                row.append(rhs[v] - lhs[v])
            row.append(rhs['cte'] - lhs['cte'])
            rows.append(row)

        elif '>=' in c:
            parts = re.split(r'>=', c)
            lhs = _parse_linear(parts[0].strip(), variables)
            rhs = _parse_linear(parts[1].strip(), variables)
            # lhs - rhs >= 0
            row = [1]
            for v in variables:
                row.append(lhs[v] - rhs[v])
            row.append(lhs['cte'] - rhs['cte'])
            rows.append(row)

        else:
            raise ValueError(f"Contrainte incomprise : '{c}'")

    return rows


# ──────────────────────────────────────────────
# Fonction principale : LBLRead
# ──────────────────────────────────────────────
def LBLRead(s):
    """
    Parse une chaîne du type "{(2i, j) | 0 <= i <= 10, i = 2j}"
    et retourne un objet LBL de pypolylib.

    Format :
      { partie_gauche | contraintes }
      - partie_gauche : expressions linéaires des variables de sortie
      - contraintes   : inégalités (>=, <=) et égalités (=) sur les variables
    """
    s = s.strip()
    if not (s.startswith('{') and s.endswith('}')):
        raise ValueError("La chaîne doit commencer par '{' et finir par '}'")
    s = s[1:-1]  # enlever { }

    # Séparer gauche | droite
    if '|' not in s:
        raise ValueError("Il manque le séparateur '|'")
    lhs_str, rhs_str = s.split('|', 1)
    lhs_str = lhs_str.strip()
    rhs_str = rhs_str.strip()

    # Détecter les variables (lettres minuscules seules ou suivies de chiffres)
    all_text = lhs_str + rhs_str
    variables = sorted(set(re.findall(r'[a-z](?![a-z])', all_text)))

    n = len(variables)  # nombre de variables (dimension de l'espace d'entrée)

    # ── Construire la matrice lattice depuis la partie gauche ──
    lhs_coeffs = _parse_lhs(lhs_str, variables)
    nb_out = len(lhs_coeffs)  # nombre de sorties

    # La matrice lattice est de taille (nb_out+1) x (n+1)
    # Dernière ligne : [0, 0, ..., 1] (ligne homogène)
    # Dernière colonne : constante
    nb_rows_lat = nb_out + 1
    nb_cols_lat = n + 1

    lat = pypolylib.MatrixAlloc(nb_rows_lat, nb_cols_lat)
    for i, coeffs in enumerate(lhs_coeffs):
        for j, v in enumerate(variables):
            pypolylib.MatrixSetValue(lat, i, j, coeffs[v])
        pypolylib.MatrixSetValue(lat, i, n, coeffs['cte'])
    # Ligne homogène
    for j in range(n):
        pypolylib.MatrixSetValue(lat, nb_out, j, 0)
    pypolylib.MatrixSetValue(lat, nb_out, n, 1)

    # ── Construire la matrice de contraintes depuis la partie droite ──
    constraint_rows = _parse_constraints(rhs_str, variables)
    nb_rows_poly = len(constraint_rows)
    nb_cols_poly = n + 2  # type + n variables + constante

    cmat = pypolylib.MatrixAlloc(nb_rows_poly, nb_cols_poly)
    for i, row in enumerate(constraint_rows):
        for j, val in enumerate(row):
            pypolylib.MatrixSetValue(cmat, i, j, val)

    poly = pypolylib.Constraints2Polyhedron(cmat, 0)

    return pypolylib.LBLAlloc(lat, poly)


# ──────────────────────────────────────────────
# Classe LBL Python
# ──────────────────────────────────────────────
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