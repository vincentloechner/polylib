import math

# 50 -> 1m points in A (~(50*2)^3) takes some time to compute...
MAX = 50

def ZP2D(lattice, i_inf, i_sup, j_inf, j_sup, k_inf, k_sup, extra="True"):
  L = set()
  for i in range(math.ceil(eval(i_inf)), math.floor(eval(i_sup))+1):
    for j in range(math.ceil(eval(j_inf)), math.floor(eval(j_sup))+1):
      for k in range(math.ceil(eval(k_inf)), math.floor(eval(k_sup))+1):
        if eval(extra):
          L.add(eval(lattice))
  return L

A = ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", "-MAX", "MAX", "-MAX", "MAX", "-MAX", "MAX")
B = ZP2D("(i, 0, k)", "0", "1000", "i/10", "i/9", "i+j/4", "i+j/3")

print("A contains", len(A), "points")
print("B contains", len(B), "points")

diff = A - B
print("diff = A-B contains", len(diff), "points")

ZP = set()


# LATTICE: 
# 4 4
#    2    0    0    0 
#    0    5    0    0 
#   52   27   66    0 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   0    0   -1   50  ]
# Inequality: [   1    0    0   50  ]
# Inequality: [  -1    0    0   50  ]
# Inequality: [   0    1    0   50  ]
# Inequality: [   0    0    1   50  ]
# Inequality: [   0   -1    0   -1  ]
ZP = ZP.union(ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", "-MAX", "MAX", "-MAX", "-1", "-MAX", "MAX"))

# Constraints 6 5
# Inequality: [   0    0   -1   50  ]
# Inequality: [   1    0    0   50  ]
# Inequality: [  -1    0    0   50  ]
# Inequality: [   0   -1    0   50  ]
# Inequality: [   0    0    1   50  ]
# Inequality: [   0    1    0   -1  ]
ZP = ZP.union(ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", "-MAX", "MAX", "1", "MAX", "-MAX", "MAX"))

# Constraints 7 5
# Inequality: [   0    0   -1   50  ]
# Inequality: [   1    0    0   50  ]
# Inequality: [  -1    0    0   50  ]
# Inequality: [   0    1    0   50  ]
# Inequality: [   0   -1    0   50  ]
# Inequality: [   0    0    1   50  ]
# Inequality: [-333 -180 -440   -1  ]
ZP = ZP.union(ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", "-MAX", "MAX", "-MAX", "MAX", "-MAX", "MAX",
                   "333*i+180*j+440*k+1<=0"))

# Constraints 7 5
# Inequality: [1348  729 1782   -1  ]
# Inequality: [   0    0   -1   50  ]
# Inequality: [   1    0    0   50  ]
# Inequality: [  -1    0    0   50  ]
# Inequality: [   0    1    0   50  ]
# Inequality: [   0   -1    0   50  ]
# Inequality: [   0    0    1   50  ]
ZP = ZP.union(ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", "-MAX", "MAX", "-MAX", "MAX", "-MAX", "MAX",
                   "1348*i+729*j+1782*k-1>=0"))

# # LATTICE:
# # 4 3
# #    2    0    0 
# #    0    0    0 
# #   52   66    0 
# #    0    0    1 
# # POLYHEDRON Dimension:2
# #            Constraints:4  Equations:0  Rays:4  Lines:0
# # Constraints 4 4
# # Inequality: [-674 -891    0  ]
# # Inequality: [   1    0   -1  ]
# # Inequality: [ 333  440    0  ]
# # Inequality: [  -1    0   36  ]
# # Rays 4 4
# # Vertex: [3960 -2997  ]/110
# # Vertex: [ 440 -333  ]/440
# # Vertex: [ 891 -674  ]/891
# # Vertex: [3564 -2696  ]/99

# ZP = ZP.union(ZP2D("(2*i, 0, 52*i+66*k)", "1", "36", "0", "0", "-333*i/440", "(-674*i)/891"))
# print("Adding holes, len(ZP) =", len(ZP))
# # UNION POLYHEDRON Dimension:2
# #            Constraints:4  Equations:0  Rays:4  Lines:0
# # Constraints 4 4
# # Inequality: [-674 -891    0  ]
# # Inequality: [  -1    0   39  ]
# # Inequality: [ 333  440    0  ]
# # Inequality: [   1    0  -38  ]
# # Rays 4 4
# # Vertex: [33858 -25615  ]/891
# # Vertex: [3861 -2921  ]/99
# # Vertex: [11583 -8762  ]/297
# # Vertex: [33858 -25612  ]/891
# ZP = ZP.union(ZP2D("(2*i, 0, 52*i+66*k)", "38", "39", "0", "0", "-333*i/440", "(-674*i)/891"))
# print("Adding holes, len(ZP) =", len(ZP))

# # UNION POLYHEDRON Dimension:2
# #            Constraints:4  Equations:0  Rays:4  Lines:0
# # Constraints 4 4
# # Inequality: [  -1    0   50  ]
# # Inequality: [-674 -891    0  ]
# # Inequality: [   1    0  -38  ]
# # Inequality: [ 674  891    3  ]
# # Rays 4 4
# # Vertex: [35640 -26963  ]/891
# # Vertex: [35640 -26960  ]/891
# # Vertex: [44550 -33700  ]/891
# # Vertex: [44550 -33703  ]/891
# ZP = ZP.union(ZP2D("(2*i, 0, 52*i+66*k)", "38", "50", "0", "0", "(-674*i-3)/891", "(-674*i)/891"))
# print("Adding holes, len(ZP) =", len(ZP))


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
print("diff2 = B-A contains", len(diff2), "points")

YP = set()
# B - A = LBL: Dimension 3

# LATTICE:
# 4 4
#    1    0    0    0 
#    0    0    0    0 
#    0    1    0    0 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [  -4    4   -1    0  ]
# Inequality: [   1    0   -9    0  ]
# Inequality: [  -1    0   10    0  ]
# Inequality: [   1    0    0 -101  ]
# Inequality: [   3   -3    1    0  ]
# Inequality: [  -1    0    0 1000  ]
#                 i    k    j
YP = YP.union(ZP2D("(i, 0, k)", "101", "1000", "i/10", "i/9", "(4*i+j)/4", "(3*i+j)/3"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0   10 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -1    0    5    0  ]
# Inequality: [ -18  -33    1  -30  ]
# Inequality: [  24   44   -1   40  ]
# Inequality: [  -1    0    0  500  ]
# Inequality: [   1    0    0  -15  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, 8*i+11*k+10)", "15", "500", "i/5", "2*i/9",
                   "(-24*i+j-40)/44", "(-18*i+j-30)/33"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    9 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -1    0    5    0  ]
# Inequality: [ -18  -33    1  -27  ]
# Inequality: [  24   44   -1   36  ]
# Inequality: [  -1    0    0  500  ]
# Inequality: [   1    0    0  -30  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, 8*i+11*k+9)", "30", "500", "i/5", "2*i/9",
                   "(-24*i+j-36)/44", "(-18*i+j-27)/33"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    8 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -1    0    5    0  ]
# Inequality: [ -18  -33    1  -24  ]
# Inequality: [  24   44   -1   32  ]
# Inequality: [  -1    0    0  500  ]
# Inequality: [   1    0    0  -19  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, 8*i+11*k+8)", "19", "500", "i/5", "2*i/9",
                   "(-24*i+j-32)/44", "(-18*i+j-24)/33"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    7 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -1    0    5    0  ]
# Inequality: [ -18  -33    1  -21  ]
# Inequality: [  24   44   -1   28  ]
# Inequality: [  -1    0    0  500  ]
# Inequality: [   1    0    0  -34  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, 8*i+11*k+7)", "34", "500", "i/5", "2*i/9",
                   "(-24*i+j-28)/44", "(-18*i+j-21)/33"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    6 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -1    0    5    0  ]
# Inequality: [ -18  -33    1  -18  ]
# Inequality: [  24   44   -1   24  ]
# Inequality: [  -1    0    0  500  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, 8*i+11*k+6)", "0", "500", "i/5", "2*i/9",
                   "(-24*i+j-24)/44", "(-18*i+j-18)/33"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    5 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -1    0    5    0  ]
# Inequality: [ -18  -33    1  -15  ]
# Inequality: [  24   44   -1   20  ]
# Inequality: [  -1    0    0  500  ]
# Inequality: [   1    0    0  -14  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, 8*i+11*k+5)", "14", "500", "i/5", "2*i/9",
                   "(-24*i+j-20)/44", "(-18*i+j-15)/33"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    4 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -1    0    5    0  ]
# Inequality: [ -18  -33    1  -12  ]
# Inequality: [  24   44   -1   16  ]
# Inequality: [  -1    0    0  500  ]
# Inequality: [   1    0    0  -29  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, 8*i+11*k+4)", "29", "500", "i/5", "2*i/9",
                   "(-24*i+j-16)/44", "(-18*i+j-12)/33"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    3 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -1    0    5    0  ]
# Inequality: [ -18  -33    1   -9  ]
# Inequality: [  24   44   -1   12  ]
# Inequality: [  -1    0    0  500  ]
# Inequality: [   1    0    0  -18  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, 8*i+11*k+3)", "18", "500", "i/5", "2*i/9",
                   "(-24*i+j-12)/44", "(-18*i+j-9)/33"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    2 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -1    0    5    0  ]
# Inequality: [ -18  -33    1   -6  ]
# Inequality: [  24   44   -1    8  ]
# Inequality: [  -1    0    0  500  ]
# Inequality: [   1    0    0  -20  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, 8*i+11*k+2)", "20", "500", "i/5", "2*i/9",
                   "(-24*i+j-8)/44", "(-18*i+j-6)/33"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    1 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   2    0   -9    0  ]  9j <= 2i
# Inequality: [  -1    0    5    0  ]  5j >= i
# Inequality: [ -18  -33    1   -3  ]  33k <= -18i+j-3
# Inequality: [  24   44   -1    4  ]  44k >= -24i+j-4
# Inequality: [  -1    0    0  500  ]
# Inequality: [   1    0    0  -35  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, 8*i+11*k+1)", "35", "500", "i/5", "2*i/9",
                   "(-24*i+j-4)/44", "(-18*i+j-3)/33"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    1    3    0    2 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   2    0   -9    0  ] 9j <= 2i
# Inequality: [  -1    0    5    0  ] 5j >= i
# Inequality: [   3   -9    1   -6  ] 9k  <= 3i + j - 6
# Inequality: [  -4   12   -1    8  ] 12k >= 4i + j - 8
# Inequality: [  -1    0    0  500  ]
# Inequality: [   1    0    0  -19  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, i+3*k+2)", "19", "500", "i/5", "2*i/9",
                   "(4*i+j-8)/12", "(3*i+j-6)/9"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    1    3    0    1 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -1    0    5    0  ]
# Inequality: [   3   -9    1   -3  ]
# Inequality: [  -4   12   -1    4  ]
# Inequality: [  -1    0    0  500  ]
# Inequality: [   1    0    0  -15  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, i+3*k+1)", "15", "500", "i/5", "2*i/9",
                   "(4*i+j-4)/12", "(3*i+j-3)/9"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    0    2    0    1 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -1    0    5    0  ]
# Inequality: [   6   -6    1   -3  ]
# Inequality: [  -8    8   -1    4  ]
# Inequality: [  -1    0    0  500  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i, 0, 2*k+1)", "0", "500", "i/5", "2*i/9",
                   "(8*i+j-4)/8", "(6*i+j-3)/6"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    1 
#    0    0    0    0 
#    0    1    0    0 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:6  Lines:0
# Constraints 5 5
# Inequality: [   2    0   -9    1  ]
# Inequality: [   6   -3    1    3  ]
# Inequality: [  -8    4   -1   -4  ]
# Inequality: [  -1    0    5   -1  ]
# Inequality: [  -1    0    0  499  ]
#                 i    k    j
YP = YP.union(ZP2D("(2*i+1, 0, k)", "0", "499", "(i+1)/5", "(2*i+1)/9",
                   "(8*i+j+4)/4", "(6*i+j+3)/3"))



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
