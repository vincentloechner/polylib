import math

def ZP2D(lattice, i_inf, i_sup, j_inf, j_sup, k_inf, k_sup, extra="True"):
  L = set()
  for i in range(math.ceil(eval(i_inf)), math.floor(eval(i_sup))+1):
    for j in range(math.ceil(eval(j_inf)), math.floor(eval(j_sup))+1):
      for k in range(math.ceil(eval(k_inf)), math.floor(eval(k_sup))+1):
        if eval(extra):
          L.add(eval(lattice))
  return L

A = ZP2D("(2*i, 0, 52*i+66*k)", "0", "39", "0", "0", "-i*333/440", "(-674*i-4)/891")
# (74, 0, 76)
# -> (37, 0, -28)
B = ZP2D("(2*i, 5*j, 52*i+27*j+66*k)", "-50", "50", "-50", "50", "-50", "50")

print("A contains", len(A), "points")
print("B contains", len(B), "points")

inter = A.intersection(B)
print("inter = A inter B contains", len(inter), "points")

ZP = set()

# LATTICE:
# 4 3
#    2    0    0
#    0    0    0
#   52   66    0
#    0    0    1
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [-674 -891   -4  ]
# Inequality: [  -1    0   39  ]
# Inequality: [ 333  440    0  ]
# Rays 3 4
# Vertex: [17160 -12987  ]/440
# Vertex: [3159 -2390  ]/81
# Vertex: [1760 -1332  ]/143
ZP = ZP.union(ZP2D("(2*i, 0, 52*i+66*k)", "0", "39", "0", "0", "-i*333/440", "(-674*i-4)/891"))
# WRONG!
# extra point == (74, 0, 76)
# from (i,j,k) = (37, 0, -28)

print("len(ZP) =", len(ZP))
print(ZP)
