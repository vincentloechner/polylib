import pypolylib_core as pl

def _lbl_repr(lbl_core):
    """
    Construit une représentation symbolique d'un LBL.
    Ex: {(3i) | 0 <= i <= 5}
    """
    parts = []
    node = lbl_core
    while node is not None:
        parts.append(_single_lbl_repr(node))
        node = node.next
    return " UNION\n".join(parts)


def _single_lbl_repr(node):
    lat = node.Lat
    poly = node.P

    nb_out = lat.nbrows - 1      # nombre de sorties (sans ligne homogène)
    nb_vars = lat.nbcolumns - 1  # nombre de variables d'entrée

    # Noms des variables : i, j, k, l, ...
    var_names = [chr(ord('i') + v) for v in range(nb_vars)]

    # ── Partie gauche : expressions de sortie ──
    out_exprs = []
    for r in range(nb_out):
        terms = []
        for c in range(nb_vars):
            coef = pl.MatrixGetValue(lat, r, c)
            if coef == 0:
                continue
            vname = var_names[c]
            if coef == 1:
                terms.append(vname)
            elif coef == -1:
                terms.append(f"-{vname}")
            else:
                terms.append(f"{coef}{vname}")
        cte = pl.MatrixGetValue(lat, r, nb_vars)
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
    if len(out_exprs) > 1 or True:
        lhs = f"({lhs})"

    # ── Partie droite : contraintes ──
    if poly is None:
        return "{" + lhs + " | <empty>}"

    nb_constraints = poly.nbconstraints
    dim = poly.dimension
    cmat = poly.constraints

    # Noms des variables du polyèdre
    poly_vars = [chr(ord('i') + v) for v in range(dim)]

    constraints = []
    for r in range(nb_constraints):
        eq_type = pl.MatrixGetValue(cmat, r, 0)  # 0=égalité, 1=inégalité
        coeffs = [pl.MatrixGetValue(cmat, r, c+1) for c in range(dim)]
        cte = pl.MatrixGetValue(cmat, r, dim+1)

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
            # Égalité : expr + cte = 0
            expr = _terms_to_str(terms)
            if cte != 0:
                constraints.append(f"{expr} = {-cte}")
            else:
                constraints.append(f"{expr} = 0")
        else:
            # Inégalité : expr + cte >= 0
            # On cherche à écrire sous forme a <= var <= b si possible
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
    if not terms:
        return "0"
    expr = terms[0]
    for t in terms[1:]:
        if t.startswith('-'):
            expr += t
        else:
            expr += "+" + t
    return expr