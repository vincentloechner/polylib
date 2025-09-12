import math


def ZP2D(lattice, i_inf, i_sup, j_inf, j_sup):
  L = set()
  for i in range(math.ceil(eval(i_inf)), math.floor(eval(i_sup))+1):
    for j in range(math.ceil(eval(j_inf)), math.floor(eval(j_sup))+1):
      # for k in range(math.ceil(eval(k_inf)), math.floor(eval(k_sup))+1):
      #   if eval(extra):
          L.add(eval(lattice))
  return L

A = ZP2D("(2*i, 2*j)", "-20", "20", "-20", "20")
B = ZP2D("(i, 0)", "0", "100", "i/4", "i/3")
print("A contains", len(A), "points")
print("B contains", len(B), "points")

diff = A - B
print("diff = A - B contains", len(diff), "points")

ZP = set()
ZP = ZP.union(ZP2D("(2, 0)", "0", "0", "0", "0"))
ZP = ZP.union(ZP2D("(2*i, 2*j)", "-20", "20", "-20", "-1"))
ZP = ZP.union(ZP2D("(2*i, 2*j)", "-20", "20", "1", "20"))
ZP = ZP.union(ZP2D("(2*i, 2*j)", "-20", "-1", "-20", "20"))

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


diff2 = B - A
print("diff2 = B - A contains", len(diff2), "points")
YP = set()
YP = YP.union(ZP2D("(2*i+1, 0)", "-20", "19", "(i+1)/2", "(2*i+1)/3"))

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
