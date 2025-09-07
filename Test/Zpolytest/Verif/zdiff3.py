import math

def ZP2D(lattice, i_inf, i_sup, j_inf, j_sup):
  L = set()
  for i in range(math.ceil(eval(i_inf)), math.floor(eval(i_sup))+1):
    for j in range(math.ceil(eval(j_inf)), math.floor(eval(j_sup))+1):
      L.add(eval(lattice))
  return L

A = ZP2D("(2*i, 2*i+j)", "1", "100", "1", "100")
B = ZP2D("(3*i, 2*i+j)", "1", "100", "1", "100")
print("A contains", len(A), "points")
print("B contains", len(B), "points")

# A - B
diff = A - B
print("diff = A-B contains", len(diff), "points")

ZP = set()
ZP = ZP.union(ZP2D("(2*i, j)", "1", "100", "(4*i+301)/3", "2*i+100"))
ZP = ZP.union(ZP2D("(2*i, j)", "1", "2", "(2*i+1)", "2*i+100"))
ZP = ZP.union(ZP2D("(200, j)", "1", "1", "201", "300"))
ZP = ZP.union(ZP2D("(6*i+2, j)", "1", "32", "6*i+3", "(4*i+101)"))
ZP = ZP.union(ZP2D("(6*i+4, j)", "1", "32", "6*i+5", "(4*i+102)"))

print("len(ZP) =", len(ZP))
if len(ZP - diff) == 0:
  print("ZP is in diff")
else:
  print("Points that should not be there:")
  print(ZP - diff)

if len(diff - ZP) == 0:
  print("diff is in ZP")
else:
  print("Missing points:")
  print(diff - ZP)

# B - A
diff2 = B - A
print("diff2 = B-A contains", len(diff2), "points")
YP = set()
YP = YP.union(ZP2D("(3*i, j)", "0", "100", "2*i+1", "3*i"))
YP = YP.union(ZP2D("(3*i, j)", "67", "100", "2*i+1", "2*i+100"))
YP = YP.union(ZP2D("(3, j)", "0", "0", "3", "102"))
YP = YP.union(ZP2D("(6*i+3, j)", "1", "32", "6*i+4", "4*i+102"))
print("len(YP) =", len(YP))
if len(YP - diff2) == 0:
  print("YP is in diff2")
else:
  print("Points that should not be there:")
  print(YP - diff2)

if len(diff2 - YP) == 0:
  print("diff2 is in YP")
else:
  print("Missing points:")
  print(diff2 - YP)