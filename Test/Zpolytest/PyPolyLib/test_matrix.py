from lbl_read import _parse_lhs, _parse_constraints, _matrix_to_polylib_string
import pypolylib, re

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
print("lat_str =")
print(lat_str)

lat = pypolylib.MatrixReadFromString(lat_str)
print("lu: nbrows=", lat.nbrows, "nbcols=", lat.nbcolumns)
for i in range(lat.nbrows):
    for j in range(lat.nbcolumns):
        print(f"  [{i}][{j}] =", pypolylib.MatrixGetValue(lat, i, j))
