def ZP2D(lattice, i_inf, i_sup, j_inf, j_sup):
  L = set()
  for i in range(eval(i_inf), eval(i_sup)+1):
    for j in range(eval(j_inf), eval(j_sup)+1):
      L.add(eval(lattice))
  return L

A = ZP2D("(2*i, 2*i+j)", "1", "100", "1", "100")
B = ZP2D("(3*i, 2*i+j)", "1", "100", "1", "100")
diff = A - B

print("A contient", len(A), "points")
print("B contient", len(B), "points")
print("diff=A-B contient", len(diff), "points")

ZP1 = ZP2D("(2,j)", "1", "1", "3", "102")
ZP2 = ZP2D("(6*i, j)", "1", "33", "((12*i+301)+2)//3", "(6*i+100)")
ZP3 = ZP2D("(6*i+4, j)", "0", "98//3", "6*i+5", "((12*i+308)+2)//3")
ZP4 = ZP2D("(6*i+2, j)", "0", "33", "6*i+3", "((12*i+304)+2)//3")

ZP = set().union(ZP1, ZP2, ZP3, ZP4)

print("len(ZP1) =", len(ZP1))
print("len(ZP2) =", len(ZP2))
print("len(ZP3) =", len(ZP3))
print("len(ZP4) =", len(ZP4))
print("len(ZP) =", len(ZP))

if ZP in diff:
  print("ZP in diff")
if diff in ZP:
  print("diff in ZP")
