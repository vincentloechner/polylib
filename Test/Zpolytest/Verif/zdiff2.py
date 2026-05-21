import math

def ZP2D(lattice, i_inf, i_sup, j_inf, j_sup):
  L = set()
  for i in range(math.ceil(eval(i_inf)), math.floor(eval(i_sup))+1):
    for j in range(math.ceil(eval(j_inf)), math.floor(eval(j_sup))+1):
      L.add(eval(lattice))
  return L

A = ZP2D("(2*i, 2*i+j)", "1", "10", "1", "10")
B = ZP2D("(2*i, 2*i+j)", "5", "15", "5", "15")

diff = A - B

print("A contient", len(A), "points")
print("B contient", len(B), "points")
print("diff=A-B contient", len(diff), "points")

# ZP1 = ZP2D("(2,j)", "1", "1", "3", "102")
# ZP2 = ZP2D("(6*i, j)", "1", "33", "((12*i+301)+2)//3", "(6*i+100)")
# ZP3 = ZP2D("(6*i+4, j)", "0", "98//3", "6*i+5", "((12*i+308)+2)//3")
# ZP4 = ZP2D("(6*i+2, j)", "0", "33", "6*i+3", "((12*i+304)+2)//3")
ZP1 = ZP2D("(2*i , j)","1","9/2","2*i+1","2*i+10")
ZP2 = ZP2D("(2*i , j)","1","10","2*i+1","2*i+4")


ZP = set().union(ZP1,ZP2)

print("len(ZP1) =", len(ZP1))
print("len(ZP2) =", len(ZP2))
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
