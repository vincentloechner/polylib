import math

# 50 -> 1m points in A (~(50*2)^3)
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
B = ZP2D("(i, 0, k)", "0", "1000", "i/10", "i/9", "j/4+i", "j/3+i")
diff = A - B

print("A contains", len(A), "points")
print("B contains", len(B), "points")
print("diff=A-B contains", len(diff), "points")


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
ZP = ZP.union(ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", "-MAX", "MAX", "-MAX", "MAX", "-MAX", "MAX", "333*i+180*j+440*k+1<=0"))

# Constraints 7 5
# Inequality: [1348  729 1782   -1  ]
# Inequality: [   1    0    0   50  ]
# Inequality: [  -1    0    0   50  ]
# Inequality: [   0    1    0   50  ]
# Inequality: [   0   -1    0   50  ]
# Inequality: [   0    0    1   50  ]
# Inequality: [   0    0   -1   50  ]
ZP = ZP.union(ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", "-MAX", "MAX", "-MAX", "MAX", "-MAX", "MAX", "1348*i+729*j+1782*k-1>=0"))


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
