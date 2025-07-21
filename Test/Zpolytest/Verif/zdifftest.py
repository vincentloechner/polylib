import math

def ZP2D(lattice, i_inf, i_sup, j_inf, j_sup):
  L = set()
  for i in range(math.ceil(eval(i_inf)), math.floor(eval(i_sup))+1):
    for j in range(math.ceil(eval(j_inf)), math.floor(eval(j_sup))+1):
      L.add(eval(lattice))
  return L

A = ZP2D("(2*i+1, j)", "-1/3", "98/3" ,"(2*i+3)/5", "(2*i+102)/5")
B = ZP2D("(2*i-1, j)",  "2/3", "101/3", "-2*i+1", "(10*i-104)/-5")


print("A contient", len(A), "points")
print("B contient", len(B), "points")


# ZP1 = ZP2D("(2,j)", "1", "1", "3", "102")
# ZP2 = ZP2D("(6*i, j)", "1", "33", "((12*i+301)+2)//3", "(6*i+100)")
# ZP3 = ZP2D("(6*i+4, j)", "0", "98//3", "6*i+5", "((12*i+308)+2)//3")
# ZP4 = ZP2D("(6*i+2, j)", "0", "33", "6*i+3", "((12*i+304)+2)//3")
# ZP1 = ZP2D("(2*i +1, j)","1","9/2","2*i+1","2*i+10")
# ZP2 = ZP2D("(2*i, j)","1","10","2*i+1","2*i+4")


ZP = set().intersection(A,B)

print("len(A) =", len(A))
print("len(B) =", len(B))
print("len(ZP) =", len(ZP))
