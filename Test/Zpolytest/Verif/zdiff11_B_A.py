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
B = ZP2D("(i, 0, k)", "0", "1000", "i/10", "i/9", "i+j/4", "i+j/3")
diff = B - A

# (34 0 35) is in B?
# i = 34, -i+10j>=0, i-9j>=0,  4k >= 4i+j, 3k <= 3i+j
# => j >= 3.4 && j <= i/9=3.77 => no integer j.
# => no

print("A contains", len(A), "points")
print("B contains", len(B), "points")
print("diff=A-B contains", len(diff), "points")


ZP = set()
i=0
# B - A = LBL: Dimension 3 
# LATTICE: 
# 4 3
#    1    0    0 
#    0    0    0 
#    0    1    0 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:4  Equations:0  Rays:4  Lines:0
# Constraints 4 4
# Inequality: [ -41   40    0  ]
# Inequality: [   1    0 -101  ]
# Inequality: [  28  -27    0  ]
# Inequality: [  -1    0 1000  ]
ZP = ZP.union(ZP2D("(i, 0, k)", "101", "1000", "0", "0", "41*i/40", "28*i/27"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    8   11   10 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [ 119  220  200  ]
# Inequality: [-160 -297 -270  ]
ZP = ZP.union(ZP2D("(2*i, 0, 8*i+11*k+10)", "0", "50", "0", "0", "(-119*i-200)/220", "(-160*i-270)/297"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    8   11    9 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [ 119  220  180  ] 220k >= -119i-180
# Inequality: [-160 -297 -243  ] 297k <= -160i-243
ZP = ZP.union(ZP2D("(2*i, 0, 8*i+11*k+9)", "0", "50", "0", "0", "(-119*i-180)/220", "(-160*i-243)/297"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# TODO:
# (34 0 35) is in ZP?
# x = 34, i=17. 0<=i<=50 ---> i=17 ok
# ---> j = 0 ok
# k >= (-119*i-180)/220=-10.01
# k <= (-160*i-243)/297=-9.97. ---> k=-10 ok
# (17, 0, -10) maps to (34, 0, 35) ok
# ---> YES but it's wrong.

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    8   11    8 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [ 119  220  160  ]
# Inequality: [-160 -297 -216  ]
ZP = ZP.union(ZP2D("(2*i, 0, 8*i+11*k+8)", "0", "50", "0", "0", "(-119*i-160)/220", "(-160*i-216)/297"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    8   11    7 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [ 119  220  140  ]
# Inequality: [-160 -297 -189  ]
ZP = ZP.union(ZP2D("(2*i, 0, 8*i+11*k+7)", "0", "50", "0", "0", "(-119*i-140)/220", "(-160*i-189)/297"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    8   11    6 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [ 119  220  120  ]
# Inequality: [-160 -297 -162  ]
ZP = ZP.union(ZP2D("(2*i, 0, 8*i+11*k+6)", "0", "50", "0", "0", "(-119*i-120)/220", "(-160*i-162)/297"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    8   11    5 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [ 119  220  100  ]
# Inequality: [-160 -297 -135  ]
ZP = ZP.union(ZP2D("(2*i, 0, 8*i+11*k+150)", "0", "50", "0", "0", "(-119*i-100)/220", "(-160*i-135)/297"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    8   11    4 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [ 119  220   80  ]
# Inequality: [-160 -297 -108  ]
ZP = ZP.union(ZP2D("(2*i, 0, 8*i+11*k+4)", "0", "50", "0", "0", "(-119*i-80)/220", "(-160*i-108)/297"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    8   11    3 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [ 119  220   60  ]
# Inequality: [-160 -297  -81  ]
ZP = ZP.union(ZP2D("(2*i, 0, 8*i+11*k+3)", "0", "50", "0", "0", "(-119*i-60)/220", "(-160*i-81)/297"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    8   11    2 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [ 119  220   40  ]
# Inequality: [-160 -297  -54  ]
ZP = ZP.union(ZP2D("(2*i, 0, 8*i+11*k+2)", "0", "50", "0", "0", "(-119*i-40)/220", "(-160*i-54)/297"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    8   11    1 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [ 119  220   20  ]
# Inequality: [-160 -297  -27  ]
ZP = ZP.union(ZP2D("(2*i, 0, 8*i+11*k+1)", "0", "50", "0", "0", "(-119*i-20)/220", "(-160*i-27)/297"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    1    3    2 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [  29  -81  -54  ]
# Inequality: [  -7   20   13  ]
ZP = ZP.union(ZP2D("(2*i, 0, i+3*k+2)", "0", "50", "0", "0", "(7*i-13)/20", "(29*i-54)/81"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    1    3    1 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [  29  -81  -27  ]
# Inequality: [  -7   20    6  ]
ZP = ZP.union(ZP2D("(2*i, 0, i+3*k+1)", "0", "50", "0", "0", "(7*i-6)/20", "(29*i-27)/81"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    0 
#    0    0    0 
#    0    2    1 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   50  ]
# Inequality: [ -41   40   20  ]
# Inequality: [  28  -27  -14  ]
ZP = ZP.union(ZP2D("(2*i, 0, 2*k+1)", "0", "50", "0", "0", "(41*i-20)/40", "(28*i-14)/27"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)

# LATTICE: 
# 4 3
#    2    0    1 
#    0    0    0 
#    0    1    0 
#    0    0    1 
# POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  56  -27   28  ]
# Inequality: [  -1    0   49  ]
# Inequality: [ -41   20  -21  ]
ZP = ZP.union(ZP2D("(2*i+1, 0, k)", "0", "49", "0", "0", "(41*i+20)/20", "(56*i-28)/27"))
print(i, len(ZP));i+=1
print("len(ZP-diff) =", len(ZP-diff), ZP-diff)



print("len(ZP) =", len(ZP))

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
