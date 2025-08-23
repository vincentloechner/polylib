import math
MAX = 100

def ZP2D(lattice, i_inf, i_sup, j_inf, j_sup, k_inf, k_sup):
  L = set()
  for i in range(math.ceil((i_inf)), math.floor((i_sup))+1):
    for j in range(math.ceil((j_inf)), math.floor((j_sup))+1):
      for k in range(math.ceil((k_inf)), math.floor((k_sup))+1):
        L.add(eval(lattice))
  return L

def zero_MAX(X):
  Y = set()
  for (x,y,z) in X:
    if x in range(0,MAX) and y in range(0,MAX) and z in range(0,MAX):
      Y.add((x,y,z))
  return Y

A = ZP2D("(i, 7*j, k)", 0, MAX, 0, MAX/7, 0, MAX) # A in (0, MAX)^3
B = ZP2D("(3*i+1, i+7*j+3, i+2*j+4*k+1)",  0, MAX/3+1, -MAX/7, MAX/7+1, -MAX/4, MAX/4+1)  # B overlaps (0, MAX)^3

A = zero_MAX(A)
B = zero_MAX(B)

AmB = A - B
BmA = B - A

print("A contient", len(A), "points")
print("B contient", len(B), "points")
print("AmB=A-B contient", len(AmB), "points")
print("BmA=B-A contient", len(BmA), "points")

# A contient 2500 points
# B contient 181 points
# AmB=A-B contient 2476 points
# BmA=B-A contient 157 points

# print("A =", A)
# print("B =", B)

# # ZP1 = ZP2D("(2,j)", "1", "1", "3", "102")
# # ZP2 = ZP2D("(6*i, j)", "1", "33", "((12*i+301)+2)//3", "(6*i+100)")
# # ZP3 = ZP2D("(6*i+4, j)", "0", "98//3", "6*i+5", "((12*i+308)+2)//3")
# # ZP4 = ZP2D("(6*i+2, j)", "0", "33", "6*i+3", "((12*i+304)+2)//3")
# ZP1 = ZP2D("(2*i +1, j)","1","9/2","2*i+1","2*i+10")
# ZP2 = ZP2D("(2*i, j)","1","10","2*i+1","2*i+4")

# ZP = set().union(ZP1,ZP2)

# print("len(ZP1) =", len(ZP1))
# print("len(ZP2) =", len(ZP2))
# print("len(ZP) =", len(ZP))

# if len(ZP - diff) == 0:
#   print("ZP is in diff")
# else:
#   print("Points that should not be there:")
#   print(ZP - diff)

# if len(diff - ZP) == 0:
#   print("diff is in ZP")
# else:
#   print("Missing points:")
#   print(diff - ZP)


# A - B:
# 5 5
#   21    0    0    0   13
#    0    7    0    0    0
#    1    2    4    0    1
#    0    0    0    1    0
#    0    0    0    0    1
ZP1 = zero_MAX(ZP2D("(21*i+13, 7*j, i+2*j+4*k+1)", 0, MAX/21+1, 0, MAX/7, -MAX, MAX/4))
# 5 5
#   21    0    0    0   13
#    0    7    0    0    0
#    1    0    2    0    0
#    0    0    0    1    0
#    0    0    0    0    1
ZP2 = zero_MAX(ZP2D("(21*i+13, 7*j, i+2*k)", 0,  MAX/21+1, 0, MAX/7, -MAX/2, MAX/2))
# 5 5
#    7    0    0    0    5
#    0    7    0    0    0
#    0    0    1    0    0
#    0    0    0    1    0
#    0    0    0    0    1
ZP3 = zero_MAX(ZP2D("(7*i+5, 7*j, k)", 0,  MAX/7+1, 0, MAX/7, 0, MAX))
# 5 5
#    7    0    0    0    4
#    0    7    0    0    0
#    0    0    1    0    0
#    0    0    0    1    0
#    0    0    0    0    1
ZP4 = zero_MAX(ZP2D("(7*i+4, 7*j, k)", 0,  MAX/7+1, 0, MAX/7, 0, MAX))
# 5 5
#    7    0    0    0    3
#    0    7    0    0    0
#    0    0    1    0    0
#    0    0    0    1    0
#    0    0    0    0    1
ZP5 = zero_MAX(ZP2D("(7*i+3, 7*j, k)", 0,  MAX/7+1, 0, MAX/7, 0, MAX))
# 5 5
#    7    0    0    0    2
#    0    7    0    0    0
#    0    0    1    0    0
#    0    0    0    1    0
#    0    0    0    0    1
ZP6 = zero_MAX(ZP2D("(7*i+2, 7*j, k)", 0,  MAX/7+1, 0, MAX/7, 0, MAX))
# 5 5
#    7    0    0    0    1
#    0    7    0    0    0
#    0    0    1    0    0
#    0    0    0    1    0
#    0    0    0    0    1
ZP7 = zero_MAX(ZP2D("(7*i+1, 7*j, k)", 0,  MAX/7+1, 0, MAX/7, 0, MAX))
# 5 5
#    7    0    0    0    0
#    0    7    0    0    0
#    0    0    1    0    0
#    0    0    0    1    0
#    0    0    0    0    1
ZP8 = zero_MAX(ZP2D("(7*i, 7*j, k)", 0,  MAX/7+1, 0, MAX/7, 0, MAX))
# 5 5
#    3    0    0    0    2
#    0    7    0    0    0
#    0    0    1    0    0
#    0    0    0    1    0
#    0    0    0    0    1
ZP9 = zero_MAX(ZP2D("(3*i+2, 7*j, k)", 0,  MAX/3+1, 0, MAX/7, 0, MAX))
# 5 5
#    3    0    0    0    0
#    0    7    0    0    0
#    0    0    1    0    0
#    0    0    0    1    0
#    0    0    0    0    1
ZP10 = zero_MAX(ZP2D("(3*i, 7*j, k)", 0,  MAX/3+1, 0, MAX/7, 0, MAX))


ZP = set().union(ZP1, ZP2, ZP3, ZP4, ZP5, ZP6, ZP7, ZP8, ZP9, ZP10)
print("len(ZP) =", len(ZP))

if AmB == ZP:
  print("AmB is equal to ZP")
else:
  print("AmB is *NOT* equal to ZP")

# ############# B - A ###################
# 5 5
#   21    0    0    0   19
#    0    7    0    0    2
#    1    2    4    0    1
#    0    0    0    1    0
#    0    0    0    0    1
YP1 = zero_MAX(ZP2D("(21*i+19, 7*j+2, i+2*j+4*k+1)", 0, MAX/21+1, 0, MAX/7+1, -MAX/4, MAX/4+1))
# 5 5
#   21    0    0    0   16
#    0    7    0    0    1
#    1    2    4    0    0
#    0    0    0    1    0
#    0    0    0    0    1
YP2 = zero_MAX(ZP2D("(21*i+16, 7*j+1, i+2*j+4*k)", 0, MAX/21+1, 0, MAX/7+1, -MAX/4, MAX/4+1))
# 5 5
#   21    0    0    0   10
#    0    7    0    0    6
#    1    2    4    0    0
#    0    0    0    1    0
#    0    0    0    0    1
YP3 = zero_MAX(ZP2D("(21*i+10, 7*j+6, i+2*j+4*k)", 0, MAX/21+1, 0, MAX/7+1, -MAX/4, MAX/4+1))
# 5 5
#   21    0    0    0    7
#    0    7    0    0    5
#    1    2    4    0    3
#    0    0    0    1    0
#    0    0    0    0    1
YP4 = zero_MAX(ZP2D("(21*i+7, 7*j+5, i+2*j+4*k+3)", 0, MAX/21+1, 0, MAX/7+1, -MAX/4, MAX/4+1))
# 5 5
#   21    0    0    0    4
#    0    7    0    0    4
#    1    2    4    0    2
#    0    0    0    1    0
#    0    0    0    0    1
YP5 = zero_MAX(ZP2D("(21*i+4, 7*j+4, i+2*j+4*k+2)", 0, MAX/21+1, 0, MAX/7+1, -MAX/4, MAX/4+1))
# 5 5
#   21    0    0    0    1
#    0    7    0    0    3
#    1    2    4    0    1
#    0    0    0    1    0
#    0    0    0    0    1
YP6 = zero_MAX(ZP2D("(21*i+1, 7*j+3, i+2*j+4*k+1)", 0, MAX/21+1, 0, MAX/7+1, -MAX/4, MAX/4+1))

YP = set().union(YP1, YP2, YP3, YP4, YP5, YP6)
print("len(YP) =", len(YP))

if BmA == YP:
  print("BmA is equal to YP")
else:
  print("BmA is *NOT* equal to YP")
