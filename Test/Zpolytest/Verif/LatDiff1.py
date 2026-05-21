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

# A (normalized) = 4 4
#    2    0    0    0
#    5   15    0    0
#    0    0    1    0
#    0    0    0    1

# B (normalized) = 4 4
#    2    0    0    0
#    0    5    0    0
#   52   27   66    0
#    0    0    0    1

A = ZP2D("(2*i, 5*i+15*j, k)", 0, MAX/2, -MAX/3-1, MAX/15+1, 0, MAX) # A in (0, MAX)^3
B = ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", 0, MAX/2, 0, MAX/5, -MAX-1, MAX/66+1)  # B overlaps (0, MAX)^3

A = zero_MAX(A)
B = zero_MAX(B)

AmB = A - B
BmA = B - A

print("A contient", len(A), "points")
print("B contient", len(B), "points")
print("AmB=A-B contient", len(AmB), "points")
print("BmA=B-A contient", len(BmA), "points")


# A - B:
# # 5 5
# #   21    0    0    0   13
# #    0    7    0    0    0
# #    1    2    4    0    1
# #    0    0    0    1    0
# #    0    0    0    0    1
# ZP1 = zero_MAX(ZP2D("(21*i+13, 7*j, i+2*j+4*k+1)", 0, MAX/21+1, 0, MAX/7, -MAX/4, MAX/4))
# 4 4
#    2    0    0    0
#    5   15    0    0
#    2    4   11   10
#    0    0    0    1
ZP = set()
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, 2*i+4*j+11*k+10)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/11+1)))
# 4 4
#    2    0    0    0
#    5   15    0    0
#    2    4   11    9
#    0    0    0    1
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, 2*i+4*j+11*k+9)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/11+1)))
# 4 4
#    2    0    0    0
#    5   15    0    0
#    2    4   11    8
#    0    0    0    1
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, 2*i+4*j+11*k+8)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/11+1)))
# 4 4
#    2    0    0    0
#    5   15    0    0
#    2    4   11    7
#    0    0    0    1
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, 2*i+4*j+11*k+7)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/11+1)))
# 4 4
#    2    0    0    0
#    5   15    0    0
#    2    4   11    6
#    0    0    0    1
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, 2*i+4*j+11*k+6)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/11+1)))
# 4 4
#    2    0    0    0
#    5   15    0    0
#    2    4   11    5
#    0    0    0    1
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, 2*i+4*j+11*k+5)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/11+1)))
# 4 4
#    2    0    0    0
#    5   15    0    0
#    2    4   11    4
#    0    0    0    1
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, 2*i+4*j+11*k+4)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/11+1)))

# 4 4
#    2    0    0    0
#    5   15    0    0
#    2    4   11    3
#    0    0    0    1
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, 2*i+4*j+11*k+3)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/11+1)))
# 4 4
#    2    0    0    0
#    5   15    0    0
#    2    4   11    2
#    0    0    0    1
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, 2*i+4*j+11*k+2)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/11+1)))
# 4 4
#    2    0    0    0
#    5   15    0    0
#    2    4   11    1
#    0    0    0    1
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, 2*i+4*j+11*k+1)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/11+1)))
# 4 4
#    2    0    0    0
#    5   15    0    0
#    1    0    3    2
#    0    0    0    1
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, i+3*k+2)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/3+1)))
# 4 4
#    2    0    0    0
#    5   15    0    0
#    1    0    3    1
#    0    0    0    1
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, i+3*k+1)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/3+1)))
# 4 4
#    2    0    0    0
#    5   15    0    0
#    1    1    2    1
#    0    0    0    1
ZP = ZP.union(zero_MAX(ZP2D("(2*i, 5*i+15*j, i+j+2*k+1)", 0, MAX/2, -MAX/3, MAX/15, -MAX, MAX/2+1)))

print("len(ZP) =", len(ZP))

if AmB == ZP:
  print("AmB is equal to ZP")
else:
  print("AmB is *NOT* equal to ZP")



# # ############# B - A ###################
# 4 4
#    2    0    0    0
#    5   15    0   10
#   13   15   66   54
#    0    0    0    1
YP1 = zero_MAX(ZP2D("(2*i, 5*i+15*j+10, 13*i+15*j+66*k+54)", 0, MAX/2, -MAX/3-1, MAX/15+1, -MAX-1, MAX/66+1))
# 4 4
#    2    0    0    0
#    5   15    0    5
#   13   15   66   27
#    0    0    0    1
YP2 = zero_MAX(ZP2D("(2*i, 5*i+15*j+5, 13*i+15*j+66*k+27)", 0, MAX/2, -MAX/3-1, MAX/15+1, -MAX-1, MAX/66+1))


# # 5 5
# #   21    0    0    0   19
# #    0    7    0    0    2
# #    1    2    4    0    1
# #    0    0    0    1    0
# #    0    0    0    0    1
# YP1 = zero_MAX(ZP2D("(21*i+19, 7*j+2, i+2*j+4*k+1)", 0, MAX/21+1, 0, MAX/7+1, -MAX/4, MAX/4+1))
# # 5 5
# #   21    0    0    0   16
# #    0    7    0    0    1
# #    1    2    4    0    0
# #    0    0    0    1    0
# #    0    0    0    0    1
# YP2 = zero_MAX(ZP2D("(21*i+16, 7*j+1, i+2*j+4*k)", 0, MAX/21+1, 0, MAX/7+1, -MAX/4, MAX/4+1))
# # 5 5
# #   21    0    0    0   10
# #    0    7    0    0    6
# #    1    2    4    0    0
# #    0    0    0    1    0
# #    0    0    0    0    1
# YP3 = zero_MAX(ZP2D("(21*i+10, 7*j+6, i+2*j+4*k)", 0, MAX/21+1, 0, MAX/7+1, -MAX/4, MAX/4+1))
# # 5 5
# #   21    0    0    0    7
# #    0    7    0    0    5
# #    1    2    4    0    3
# #    0    0    0    1    0
# #    0    0    0    0    1
# YP4 = zero_MAX(ZP2D("(21*i+7, 7*j+5, i+2*j+4*k+3)", 0, MAX/21+1, 0, MAX/7+1, -MAX/4, MAX/4+1))
# # 5 5
# #   21    0    0    0    4
# #    0    7    0    0    4
# #    1    2    4    0    2
# #    0    0    0    1    0
# #    0    0    0    0    1
# YP5 = zero_MAX(ZP2D("(21*i+4, 7*j+4, i+2*j+4*k+2)", 0, MAX/21+1, 0, MAX/7+1, -MAX/4, MAX/4+1))
# # 5 5
# #   21    0    0    0    1
# #    0    7    0    0    3
# #    1    2    4    0    1
# #    0    0    0    1    0
# #    0    0    0    0    1
# YP6 = zero_MAX(ZP2D("(21*i+1, 7*j+3, i+2*j+4*k+1)", 0, MAX/21+1, 0, MAX/7+1, -MAX/4, MAX/4+1))

# YP = set().union(YP1, YP2, YP3, YP4, YP5, YP6)
YP = set().union(YP1, YP2)
print("len(YP) =", len(YP))

if BmA == YP:
  print("BmA is equal to YP")
else:
  print("BmA is *NOT* equal to YP")
  print("BmA - YP =", BmA - YP)
  print("YP - BmA =", YP - BmA)
