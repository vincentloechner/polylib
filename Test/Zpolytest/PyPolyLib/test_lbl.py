import pypolylib
from lbl_read import LBLRead, _parse_lhs, _parse_constraints, _matrix_to_polylib_string
import re

s = "(i+j) | 0 <= i <= 10, i = 2j"
lhs_str, rhs_str = s.split('|', 1)
variables = sorted(set(re.findall(r'[a-z](?![a-z])', lhs_str + rhs_str)))
n = len(variables)

lhs_coeffs = _parse_lhs(lhs_str.strip(), variables)
nb_out = len(lhs_coeffs)
nb_rows_lat = nb_out + 1
nb_cols_lat = n + 1

lat_values = []
for coeffs in lhs_coeffs:
    for v in variables:
        lat_values.append(coeffs[v])
    lat_values.append(coeffs['cte'])
for j in range(n):
    lat_values.append(0)
lat_values.append(1)
lat_str = _matrix_to_polylib_string(nb_rows_lat, nb_cols_lat, lat_values)
lat = pypolylib.MatrixReadFromString(lat_str)

constraint_rows = _parse_constraints(rhs_str.strip(), variables)
nb_rows_poly = len(constraint_rows)
nb_cols_poly = n + 2
poly_values = []
for row in constraint_rows:
    poly_values.extend(row)
cmat_str = _matrix_to_polylib_string(nb_rows_poly, nb_cols_poly, poly_values)
cmat = pypolylib.MatrixReadFromString(cmat_str)

print("cmat_str =")
print(cmat_str)
print("cmat lu: nbrows=", cmat.nbrows, "nbcols=", cmat.nbcolumns)

poly = pypolylib.Constraints2Polyhedron(cmat, 0)
print("poly dimension=", poly.dimension)

a = pypolylib.LBLAlloc(lat, poly)
print("LBL cree")

lat2 = a.Lat
print("Lat apres LBLAlloc: nbrows=", lat2.nbrows, "nbcols=", lat2.nbcolumns)
for i in range(lat2.nbrows):
    for j in range(lat2.nbcolumns):
        print(f"  [{i}][{j}] =", pypolylib.MatrixGetValue(lat2, i, j))
