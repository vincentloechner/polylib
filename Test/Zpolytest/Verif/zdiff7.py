import math

def ZP2D(lattice, i_inf, i_sup, j_inf, j_sup):
  L = set()
  for i in range(math.ceil(eval(i_inf)), math.floor(eval(i_sup))+1):
    for j in range(math.ceil(eval(j_inf)), math.floor(eval(j_sup))+1):
      L.add(eval(lattice))
  return L

A = ZP2D("(2*i, j)", "1", "10", "1", "10")
B = ZP2D("(3*i, j)", "1", "20", "5", "5")
print("A contains", len(A), "points")
print("B contains", len(B), "points")

# A - B
diff = A - B
print("diff = A-B contains", len(diff), "points")

ZP = set()
ZP = ZP.union(ZP2D("(2*i, j)", "1", "10", "1", "4"))
ZP = ZP.union(ZP2D("(2*i, j)", "1", "10", "6", "10"))
ZP = ZP.union(ZP2D("(2*i, j)", "1", "2", "1", "10"))
ZP = ZP.union(ZP2D("(20, j)", "1", "1", "1", "10"))
ZP = ZP.union(ZP2D("(6*i+2, 5)", "0", "3", "1", "1"))
ZP = ZP.union(ZP2D("(6*i+4, 5)", "0", "2", "1", "1"))

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
YP = YP.union(ZP2D("(3*i, 5)", "7", "20", "1", "1"))
YP = YP.union(ZP2D("(3, 5)", "1", "1", "1", "1"))
YP = YP.union(ZP2D("(6*i+3, 5)", "0", "9", "1", "1"))

print("len(ZP) =", len(YP))
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
