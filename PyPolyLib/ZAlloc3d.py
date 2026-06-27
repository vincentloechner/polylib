from pypolylib import LBLRead

# Test inspiré de ZAlloc3d.in (exemple IMPACT paper)
print("=== ZAlloc3d ===")

a = LBLRead("{(i+j) | 0 <= i <= 10, i = 2j}")
print("a =", a)
# attendu : {(3i) | 0 <= i <= 5}

expected = LBLRead("{(3i) | 0 <= i <= 5}")
check1 = a.included(expected)
check2 = expected.included(a)
print("Check (a == expected):", check1 and check2)
assert check1 and check2, "ERREUR : résultat inattendu !"
print("OK")
