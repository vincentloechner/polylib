"""
pypolylib.py — Interface Python pour la bibliothèque PolyLib (Z-polyèdres).

Ce module expose les types LBL (Lattice-Based Lattice) et Transfo
pour manipuler des ensembles de points entiers (Z-polyèdres) en Python.

Utilisation typique :
    from pypolylib import LBLRead, Transfo
    a = LBLRead("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
    f = Transfo("(i,j -> 2i+1, i+3j)")
    print(a.image(f))
"""

import pypolylib_core as pl
import re

from lbl_repr import _lbl_repr


# ──────────────────────────────────────────────────────────────
# Fonctions internes de parsing
# ──────────────────────────────────────────────────────────────

def _parse_linear(expr, variables):
    """
    Parse une expression linéaire en termes des variables données.

    Args:
        expr (str): Expression linéaire, ex: "2i+3j-1"
        variables (list): Liste des noms de variables, ex: ['i', 'j']

    Returns:
        dict: Coefficients par variable + constante, ex: {'i':2, 'j':3, 'cte':-1}
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
                raise ValueError(f"Terme incompris : '{term}' dans '{expr}'")
    return coeffs


def _parse_lhs(lhs_str, variables):
    """
    Parse la partie gauche d'un LBL : "(expr1, expr2, ...)".

    Args:
        lhs_str (str): Partie gauche, ex: "(i+j)" ou "(2i, j+1)"
        variables (list): Liste des noms de variables

    Returns:
        list: Liste de dicts de coefficients, un par expression de sortie
    """
    lhs_str = lhs_str.strip()
    if lhs_str.startswith('(') and lhs_str.endswith(')'):
        lhs_str = lhs_str[1:-1]
    exprs = [e.strip() for e in lhs_str.split(',')]
    return [_parse_linear(e, variables) for e in exprs]


def _parse_constraints(rhs_str, variables):
    """
    Parse les contraintes d'un LBL en lignes au format PolyLib.

    Supporte : >=, <=, =, >, <, doubles inégalités (ex: 0 <= i <= 10).
    Chaque ligne retournée a la forme [type, coeff_v1, ..., coeff_vn, constante]
    où type=0 pour égalité et type=1 pour inégalité (>= 0).

    Args:
        rhs_str (str): Partie droite, ex: "0 <= i <= 10, i = 2j"
        variables (list): Liste des noms de variables

    Returns:
        list: Liste de listes d'entiers représentant les contraintes PolyLib
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
            raise ValueError(f"Contrainte incomprise : '{c}'")
    return rows


def _matrix_to_polylib_string(nb_rows, nb_cols, values):
    """
    Construit une chaîne au format PolyLib pour Matrix_Read.

    Args:
        nb_rows (int): Nombre de lignes
        nb_cols (int): Nombre de colonnes
        values (list): Valeurs de la matrice, ligne par ligne

    Returns:
        str: Chaîne au format PolyLib, ex: "2 3\\n1 0 0\\n0 1 0\\n"
    """
    s = f"{nb_rows} {nb_cols}\n"
    idx = 0
    for i in range(nb_rows):
        s += " ".join(str(values[idx + j]) for j in range(nb_cols)) + "\n"
        idx += nb_cols
    return s


def LBLRead(s):
    """
    Crée un LBL depuis une chaîne symbolique.

    C'est un raccourci pour LBL(s).

    Args:
        s (str): Chaîne au format "{(expr1, ...) | contraintes}"

    Returns:
        LBL: L'objet LBL correspondant

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
    Z-polyèdre (Lattice-Based Lattice) — ensemble de points entiers.

    Un LBL est défini par un lattice (matrice de transformation) et
    un polyèdre de contraintes. Il représente l'image d'un ensemble
    de points entiers par une fonction linéaire affine.

    Création :
        a = LBL("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")
        # ou :
        a = LBLRead("{(i, j) | 1 <= i <= 10, 1 <= j <= 10}")

    Opérations :
        a + b   # union
        a - b   # différence
        a * b   # intersection
    """

    def _LBLRead(self, s):
        """
        Parse une chaîne symbolique et construit le LBL interne.

        Args:
            s (str): Chaîne au format "{(expr1, ...) | contraintes}"
                     Les variables sont des lettres minuscules (i, j, k, ...).
                     Contraintes supportées : <=, >=, =, <, >, doubles inégalités.

        Returns:
            pypolylib_core.LBL: L'objet LBL C correspondant
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

        # Matrice lattice
        lhs_coeffs = _parse_lhs(lhs_str, variables)
        nb_out = len(lhs_coeffs)
        nb_rows_lat = nb_out + 1
        nb_cols_lat = n + 1

        lat_values = []
        for coeffs in lhs_coeffs:
            for v in variables:
                lat_values.append(coeffs[v])
            lat_values.append(coeffs['cte'])
        # ligne homogène
        for j in range(n):
            lat_values.append(0)
        lat_values.append(1)

        lat_str = _matrix_to_polylib_string(nb_rows_lat, nb_cols_lat, lat_values)
        lat = pl.MatrixReadFromString(lat_str)

        # Matrice de contraintes
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
        Initialise un LBL, optionnellement depuis une chaîne symbolique.

        Args:
            s (str, optional): Chaîne au format "{(expr) | contraintes}".
                               Si None, crée un LBL vide.
        """
        self._lbl = None
        if s is not None:
            self._lbl = self._LBLRead(s)

    def Print(self):
        """Affiche le LBL sous forme symbolique."""
        print(self.__repr__())

    def __repr__(self):
        """Retourne une représentation symbolique du LBL, ex: {(3i) | i >= 0, i <= 5}"""
        if self._lbl is None:
            return "LBL(vide)"
        return _lbl_repr(self._lbl)

    def intersection(self, other):
        """
        Calcule l'intersection de ce LBL avec un autre.

        Args:
            other (LBL): L'autre LBL (même dimension)

        Returns:
            LBL: L'intersection des deux LBL
        """
        result = LBL()
        result._lbl = self._lbl.intersection(other._lbl)
        return result

    def difference(self, other):
        """
        Calcule la différence de ce LBL avec un autre (self - other).

        Args:
            other (LBL): L'autre LBL (même dimension)

        Returns:
            LBL: La différence des deux LBL
        """
        result = LBL()
        result._lbl = self._lbl.difference(other._lbl)
        return result

    def union(self, other):
        """
        Calcule l'union de ce LBL avec un autre.

        Args:
            other (LBL): L'autre LBL (même dimension)

        Returns:
            LBL: L'union des deux LBL
        """
        result = LBL()
        result._lbl = self._lbl.union(other._lbl)
        return result

    def included(self, other):
        """
        Teste si ce LBL est inclus dans un autre.

        Args:
            other (LBL): L'autre LBL

        Returns:
            bool: True si self ⊆ other, False sinon
        """
        return self._lbl.included(other._lbl)

    def zdomain(self):
        """
        Calcule le domaine entier (Z-domain) de ce LBL,
        en éliminant les variables existentielles.

        Returns:
            LBL: Le Z-domain correspondant
        """
        result = LBL()
        result._lbl = self._lbl.zdomain()
        return result

    def image(self, transfo):
        """
        Calcule l'image de ce LBL par une transformation affine.

        Args:
            transfo (Transfo): La transformation à appliquer

        Returns:
            LBL: L'image du LBL par la transformation

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
        Calcule la préimage de ce LBL par une transformation affine.

        Args:
            transfo (Transfo): La transformation à inverser

        Returns:
            LBL: La préimage du LBL par la transformation
        """
        result = LBL()
        result._lbl = pl.LBLPreimage(self._lbl, transfo._mat)
        return result

    def __add__(self, other):
        """Union : a + b  ≡  a.union(b)"""
        return self.union(other)

    def __sub__(self, other):
        """Différence : a - b  ≡  a.difference(b)"""
        return self.difference(other)

    def __mul__(self, other):
        """Intersection : a * b  ≡  a.intersection(b)"""
        return self.intersection(other)
    
    def plot(self):
        """Affiche le LBL avec matplotlib (dimension 2 uniquement)."""
        from lbl_plot import lbl_plot
        lbl_plot(self)

    def __iter__(self):
        """
        Énumère tous les points entiers du LBL.
        Parcourt les points entiers du polyèdre et calcule leur image par le lattice.
        Lève une erreur si le polyèdre n'est pas borné.

        Yields:
            tuple: Coordonnées du point entier dans l'image du lattice
        """
        seen = set()
        node = self._lbl

        while node is not None:
            lat = node.Lat
            poly = node.P

            if poly is None:
                node = node.next
                continue

            nb_out = lat.nbrows - 1
            nb_vars = lat.nbcolumns - 1
            dim = poly.dimension
            nb_constraints = poly.nbconstraints
            cmat = poly.constraints

            # Vérifier si le polyèdre est borné : NbBid == 0 et pas de ray infini
            if poly.nbbid > 0:
                raise ValueError("Le polyèdre n'est pas borné, énumération impossible")

            # Extraire les contraintes : coeffs * x + cte >= 0
            constraints = []
            for r in range(nb_constraints):
                eq_type = pl.MatrixGetValue(cmat, r, 0)
                coeffs = [pl.MatrixGetValue(cmat, r, c+1) for c in range(dim)]
                cte = pl.MatrixGetValue(cmat, r, dim+1)
                constraints.append((eq_type, coeffs, cte))

            # Trouver les bornes depuis les rayons/vertices
            bounds = []
            for v in range(dim):
                vals = []
                for r in range(poly.nbrays):
                    # accès aux rayons via les contraintes min/max
                    pass
                # Méthode alternative : chercher bornes depuis les contraintes
                lo, hi = -10000, 10000
                for (eq_type, coeffs, cte) in constraints:
                    if coeffs[v] != 0 and all(coeffs[w] == 0 for w in range(dim) if w != v):
                        if coeffs[v] > 0:
                            lo = max(lo, -cte // coeffs[v])
                        else:
                            hi = min(hi, cte // (-coeffs[v]))
                bounds.append((lo, hi))

            # Générer tous les points entiers dans la bounding box
            def _enumerate(idx, point):
                if idx == dim:
                    # Vérifier toutes les contraintes
                    for (eq_type, coeffs, cte) in constraints:
                        val = sum(coeffs[c] * point[c] for c in range(dim)) + cte
                        if val < 0:
                            return
                        if eq_type == 0 and val != 0:
                            return
                    # Calculer l'image par le lattice
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
                lo, hi = bounds[idx]
                for val in range(lo, hi+1):
                    yield from _enumerate(idx+1, point + [val])

            yield from _enumerate(0, [])
            node = node.next


# ──────────────────────────────────────────────────────────────
# Classe Transfo
# ──────────────────────────────────────────────────────────────

class Transfo:
    """
    Transformation linéaire affine sur des points entiers.

    Une Transfo est définie par une notation symbolique de la forme
    "(i,j -> 3i+1, 2i+5j)", qui est convertie en matrice PolyLib.

    Utilisation :
        f = Transfo("(i,j -> 3i+1, 2i+5j)")
        b = a.image(f)      # image d'un LBL par f
        c = a.preimage(f)   # préimage d'un LBL par f
    """

    def __init__(self, s=None):
        """
        Initialise une Transfo depuis une chaîne symbolique.

        Args:
            s (str, optional): Chaîne au format "(var1,var2,... -> expr1, expr2, ...)".
                               Si None, crée une Transfo vide.
        """
        self._mat = None
        if s is not None:
            self._mat = self._parse(s)

    def _parse(self, s):
        """
        Parse la chaîne symbolique et construit la matrice interne.

        Args:
            s (str): Chaîne au format "(i,j -> 3i+1, 2i+5j)"

        Returns:
            pypolylib_core.Matrix: La matrice de transformation PolyLib
        """
        s = s.strip()
        if s.startswith('(') and s.endswith(')'):
            s = s[1:-1]
        if '->' not in s:
            raise ValueError("Il manque '->' dans la transformation")

        lhs_str, rhs_str = s.split('->', 1)

        # Variables d'entrée : i, j, k, ...
        variables = [v.strip() for v in lhs_str.split(',')]
        n = len(variables)

        # Expressions de sortie
        out_exprs = [e.strip() for e in rhs_str.split(',')]
        nb_out = len(out_exprs)

        # Matrice (nb_out+1) x (n+1)
        nb_rows = nb_out + 1
        nb_cols = n + 1

        values = []
        for expr in out_exprs:
            coeffs = _parse_linear(expr, variables)
            for v in variables:
                values.append(coeffs[v])
            values.append(coeffs['cte'])
        # ligne homogène
        for j in range(n):
            values.append(0)
        values.append(1)

        mat_str = _matrix_to_polylib_string(nb_rows, nb_cols, values)
        return pl.MatrixReadFromString(mat_str)

    def __repr__(self):
        """
        Retourne la représentation symbolique de la transformation.

        Returns:
            str: Ex: "(i, j -> 3i+1, 2i+5j)"
        """
        if self._mat is None:
            return "Transfo(vide)"
        nb_out = self._mat.nbrows - 1
        nb_vars = self._mat.nbcolumns - 1
        var_names = [chr(ord('i') + v) for v in range(nb_vars)]

        lhs = ", ".join(var_names)

        from lbl_repr import _terms_to_str
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
        Compose deux transformations : self * other.

        Args:
            other (Transfo): L'autre transformation

        Returns:
            Transfo: La composition des deux transformations
        """
        result = Transfo()
        result._mat = pl.MatrixProduct(self._mat, other._mat)
        return result

    def inverse(self):
        """
        Calcule l'inverse de cette transformation.

        Returns:
            Transfo: La transformation inverse

        Raises:
            RuntimeError: Si la matrice n'est pas inversible
        """
        result = Transfo()
        result._mat = pl.MatrixInverse(self._mat)
        return result

    def __mul__(self, other):
        """Composition : f * g  ≡  f.compose(g)"""
        return self.compose(other)