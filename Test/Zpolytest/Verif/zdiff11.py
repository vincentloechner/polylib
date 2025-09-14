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
# 4 3
#    2    0    0 
#    0    0    0 
#   52   66    0 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 4
# Inequality: [-674 -891    0  ]
# Inequality: [ 674  891    3  ]
# Inequality: [  -1    0   39  ]
# Inequality: [ 333  440    0  ]
# Inequality: [   1    0   -1  ] (BInter) or [  26   33   -1  ] (B)
# Rays 5 4
# Vertex: [ 440 -333  ]/440
# Vertex: [3861 -2921  ]/99
# Vertex: [1320 -999  ]/143
# Vertex: [11583 -8762  ]/297
# Vertex: [ 891 -674  ]/891
ZP = ZP.union(ZP2D("(2*i, 0, 52*i+66*k)", "1", "39", "0", "0", "(-674*i-3)/891", "(-674*i)/891",
                   "333*i+440*k>=0"))
print("len(ZP) =", len(ZP))
if len(ZP - diff) == 0:
  print("ok")

# UNION POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [-674 -891   -4  ]
# Inequality: [ 333  440    0  ]
# Inequality: [  -1    0   39  ]
# Rays 3 4
# Vertex: [17160 -12987  ]/440
# Vertex: [1760 -1332  ]/143
# Vertex: [3159 -2390  ]/81
ZP = ZP.union(ZP2D("(2*i, 0, 52*i+66*k)", "-50", "39", "0", "0", "-i*333/440", "(-674*i-4)/891"))
# WRONG!
# extra point == (74, 0, 76)
# from (i,j,k) = (37, 0, -28)

print("len(ZP) =", len(ZP))
if len(ZP - diff) == 0:
  print("ZP is in diff")
else:
  print("Points that should not be there:")
  print(ZP - diff)
# UNION POLYHEDRON Dimension:2
#            Constraints:4  Equations:0  Rays:4  Lines:0
# Constraints 4 4
# Inequality: [  -1    0   50  ]
# Inequality: [-674 -891    0  ]
# Inequality: [   1    0  -40  ]
# Inequality: [ 674  891    3  ]
# Rays 4 4
# Vertex: [35640 -26963  ]/891
# Vertex: [35640 -26960  ]/891
# Vertex: [44550 -33700  ]/891
# Vertex: [44550 -33703  ]/891
ZP = ZP.union(ZP2D("(2*i, 0, 52*i+66*k)", "40", "50", "0", "0", "(-674*i-3)/891", "(-674*i)/891"))

# LATTICE: 
# 4 4
#    2    0    0    0 
#    0    5    0    0 
#   52   27   66    0 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [   0   -1    0   -1  ]
# Inequality: [   1    0    0   50  ]
# Inequality: [  -1    0    0   50  ]
# Inequality: [   0    1    0   50  ]
# Inequality: [   0    0    1   50  ]
# Inequality: [   0    0   -1   50  ]
ZP = ZP.union(ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", "-MAX", "MAX", "-MAX", "-1", "-MAX", "MAX"))

# Constraints 6 5
# Inequality: [   0    1    0   -1  ]
# Inequality: [   1    0    0   50  ]
# Inequality: [  -1    0    0   50  ]
# Inequality: [   0   -1    0   50  ]
# Inequality: [   0    0    1   50  ]
# Inequality: [   0    0   -1   50  ]
ZP = ZP.union(ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", "-MAX", "MAX", "1", "MAX", "-MAX", "MAX"))

# Constraints 7 5
# Inequality: [-333 -180 -440   -1  ]
# Inequality: [   1    0    0   50  ]
# Inequality: [  -1    0    0   50  ]
# Inequality: [   0    1    0   50  ]
# Inequality: [   0   -1    0   50  ]
# Inequality: [   0    0    1   50  ]
# Inequality: [   0    0   -1   50  ]
ZP = ZP.union(ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", "-MAX", "MAX", "-MAX", "MAX", "-MAX", "MAX",
                   "333*i+180*j+440*k+1<=0"))

# Constraints 7 5
# Inequality: [1348  729 1782   -1  ]
# Inequality: [   1    0    0   50  ]
# Inequality: [  -1    0    0   50  ]
# Inequality: [   0    1    0   50  ]
# Inequality: [   0   -1    0   50  ]
# Inequality: [   0    0    1   50  ]
# Inequality: [   0    0   -1   50  ]
ZP = ZP.union(ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", "-MAX", "MAX", "-MAX", "MAX", "-MAX", "MAX",
                   "1348*i+729*j+1782*k-1>=0"))


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
# YP = YP.union(ZP2D("(2*i, 0, 52*i+66*k)", "1", "39", "0", "0", "(-674*i-3)/891", "(-674*i)/891",                   "333*i+440*k>=0"))
