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
diff = A - B

print("A contains", len(A), "points")
print("B contains", len(B), "points")
print("diff=A-B contains", len(diff), "points")

print("diff=B-A contains", len(B-A), "points")

ZP = set()
ZP = ZP.union(ZP2D("(2, 0)", "0", "0", "0", "0"))
ZP = ZP.union(ZP2D("(2*i, 2*j)", "-20", "20", "-20", "-1"))
ZP = ZP.union(ZP2D("(2*i, 2*j)", "-20", "20", "1", "20"))
if len(ZP - diff) == 0:
  print("ok until here")
ZP = ZP.union(ZP2D("(2*i, 2*j)", "11", "20", "-20", "20"))
if len(ZP - diff) == 0:
  print("ZP is in diff")
else:
  print("problem")
  print(ZP - diff)
ZP = ZP.union(ZP2D("(2*i, 2*j)", "-20", "-1", "-20", "20"))
if len(ZP - diff) == 0:
  print("ZP is in diff")
else:
  print("Points that should not be there:")
  print(ZP - diff)



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
