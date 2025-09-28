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
# Inequality: [ 333  440    0  ]
# Inequality: [   1    0   -1  ]
# Inequality: [ 674  891    3  ]
# Inequality: [  -1    0   36  ]
# Rays 5 4
# Vertex: [10692 -8089  ]/297
# Vertex: [ 440 -333  ]/440
# Vertex: [1320 -999  ]/143
# Vertex: [ 891 -674  ]/891
# Vertex: [3564 -2696  ]/99
ZP = ZP.union(ZP2D("(2*i, 0, 52*i+66*k)", "1", "36", "0", "0", "(-674*i-3)/891", "(-674*i)/891",
                   "333*i+440*k>=0"))
print("Adding holes, len(ZP) =", len(ZP))
# UNION POLYHEDRON Dimension:2
#            Constraints:4  Equations:0  Rays:4  Lines:0
# Constraints 4 4
# Inequality: [-674 -891    0  ]
# Inequality: [  -1    0   39  ]
# Inequality: [ 674  891    3  ]
# Inequality: [   1    0  -38  ]
# Rays 4 4
# Vertex: [33858 -25615  ]/891
# Vertex: [3861 -2921  ]/99
# Vertex: [11583 -8762  ]/297
# Vertex: [33858 -25612  ]/891
ZP = ZP.union(ZP2D("(2*i, 0, 52*i+66*k)", "38", "39", "0", "0", "(-674*i-3)/891", "(-674*i)/891"))
print("Adding holes, len(ZP) =", len(ZP))
# UNION POLYHEDRON Dimension:2
#            Constraints:3  Equations:0  Rays:3  Lines:0
# Constraints 3 4
# Inequality: [  -1    0   36  ]
# Inequality: [ 333  440    0  ]
# Inequality: [-674 -891   -4  ]
# Rays 3 4
# Vertex: [1760 -1332  ]/143
# Vertex: [3960 -2997  ]/110
# Vertex: [32076 -24268  ]/891
ZP = ZP.union(ZP2D("(2*i, 0, 52*i+66*k)", "-50", "36", "0", "0", "-i*333/440", "(-674*i-4)/891"))
print("Adding holes, len(ZP) =", len(ZP))
# UNION POLYHEDRON Dimension:2
#            Constraints:4  Equations:0  Rays:4  Lines:0
# Constraints 4 4
# Inequality: [ 333  440    0  ]
# Inequality: [  -1    0   39  ]
# Inequality: [-674 -891   -4  ]
# Inequality: [   1    0  -38  ]
# Rays 4 4
# Vertex: [33858 -25616  ]/891
# Vertex: [3159 -2390  ]/81
# Vertex: [17160 -12987  ]/440
# Vertex: [8360 -6327  ]/220
ZP = ZP.union(ZP2D("(2*i, 0, 52*i+66*k)", "38", "39", "0", "0", "-i*333/440", "(-674*i-4)/891"))
print("Adding holes, len(ZP) =", len(ZP))
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
print("Adding holes, len(ZP) =", len(ZP))

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


# LATTICE:
# 4 4
#    1    0    0    0 
#    0    0    0    0 
#    0    1    0    0 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [  -4    4   -1    0  ] depends on k...
# Inequality: [   1    0   -9    0  ] 
# Inequality: [  -1    0   10    0  ] 
# Inequality: [   1    0    0 -101  ] x
# Inequality: [   3   -3    1    0  ]
# Inequality: [  -1    0    0 1000  ] x
YP = YP.union(ZP2D("(i, 0, k)", "101", "1000", "0", "0", "max(150*i+198*j, i/5)", "min(200*i+264*j, 2*i/9)"))

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0   10 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [  -1    0    0   50  ]
# Inequality: [ -18  -33    1  -30  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  24   44   -1   40  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0  -10    0  ]/11
# Vertex: [1650 -920  330  ]/33
# Vertex: [14850 -8270 3300  ]/297
# Vertex: [4950 -2765 1100  ]/99
# Vertex: [1100 -615  220  ]/22

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    9 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [  -1    0    0   50  ]
# Inequality: [ -18  -33    1  -27  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  24   44   -1   36  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0   -9    0  ]/11
# Vertex: [1650 -917  330  ]/33
# Vertex: [14850 -8243 3300  ]/297
# Vertex: [4950 -2756 1100  ]/99
# Vertex: [1100 -613  220  ]/22

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    8 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [  -1    0    0   50  ]
# Inequality: [ -18  -33    1  -24  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  24   44   -1   32  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0   -8    0  ]/11
# Vertex: [1650 -914  330  ]/33
# Vertex: [14850 -8216 3300  ]/297
# Vertex: [4950 -2747 1100  ]/99
# Vertex: [1100 -611  220  ]/22

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    7 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [  -1    0    0   50  ]
# Inequality: [ -18  -33    1  -21  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  24   44   -1   28  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0   -7    0  ]/11
# Vertex: [1650 -911  330  ]/33
# Vertex: [14850 -8189 3300  ]/297
# Vertex: [4950 -2738 1100  ]/99
# Vertex: [1100 -609  220  ]/22

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
# Inequality: [  -1    0    0   50  ]
# Inequality: [ -18  -33    1  -18  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  24   44   -1   24  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0   -6    0  ]/11
# Vertex: [1650 -908  330  ]/33
# Vertex: [1350 -742  300  ]/27
# Vertex: [4950 -2729 1100  ]/99
# Vertex: [1100 -607  220  ]/22

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    5 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [  -1    0    0   50  ]
# Inequality: [ -18  -33    1  -15  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  24   44   -1   20  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0   -5    0  ]/11
# Vertex: [1650 -905  330  ]/33
# Vertex: [14850 -8135 3300  ]/297
# Vertex: [4950 -2720 1100  ]/99
# Vertex: [ 100  -55   20  ]/2

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    4 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [  -1    0    0   50  ]
# Inequality: [ -18  -33    1  -12  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  24   44   -1   16  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0   -4    0  ]/11
# Vertex: [ 150  -82   30  ]/3
# Vertex: [14850 -8108 3300  ]/297
# Vertex: [4950 -2711 1100  ]/99
# Vertex: [1100 -603  220  ]/22

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    3 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [  -1    0    0   50  ]
# Inequality: [ -18  -33    1   -9  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  24   44   -1   12  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0   -3    0  ]/11
# Vertex: [1650 -899  330  ]/33
# Vertex: [14850 -8081 3300  ]/297
# Vertex: [4950 -2702 1100  ]/99
# Vertex: [1100 -601  220  ]/22

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    2 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [  -1    0    0   50  ]
# Inequality: [ -18  -33    1   -6  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  24   44   -1    8  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0   -2    0  ]/11
# Vertex: [1650 -896  330  ]/33
# Vertex: [14850 -8054 3300  ]/297
# Vertex: [4950 -2693 1100  ]/99
# Vertex: [1100 -599  220  ]/22

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    8   11    0    1 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [  -1    0    0   50  ]
# Inequality: [ -18  -33    1   -3  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  24   44   -1    4  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0   -1    0  ]/11
# Vertex: [1650 -893  330  ]/33
# Vertex: [14850 -8027 3300  ]/297
# Vertex: [ 450 -244  100  ]/9
# Vertex: [1100 -597  220  ]/22

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    1    3    0    2 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [  -1    0    0   50  ]
# Inequality: [   3   -9    1   -6  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -4   12   -1    8  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0   -2    0  ]/3
# Vertex: [ 450  154   90  ]/9
# Vertex: [4050 1396  900  ]/81
# Vertex: [1350  457  300  ]/27
# Vertex: [ 300  101   60  ]/6

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#    1    3    0    1 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:5  Lines:0
# Constraints 5 5
# Inequality: [  -1    0    0   50  ]
# Inequality: [   3   -9    1   -3  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -4   12   -1    4  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0   -1    0  ]/3
# Vertex: [ 450  157   90  ]/9
# Vertex: [4050 1423  900  ]/81
# Vertex: [1350  466  300  ]/27
# Vertex: [ 300  103   60  ]/6

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
# Inequality: [  -1    0    0   50  ]
# Inequality: [   6   -6    1   -3  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -8    8   -1    4  ]
# Inequality: [  -1    0    5    0  ]
# Rays 5 5
# Vertex: [   0   -1    0  ]/2
# Vertex: [ 300  307   60  ]/6
# Vertex: [2700 2773  600  ]/54
# Vertex: [ 450  458  100  ]/9
# Vertex: [ 200  203   40  ]/4

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
# Inequality: [  -1    0    0   49  ]
# Inequality: [  -1    0    5   -1  ]
# Inequality: [  -8    4   -1   -4  ]
# Inequality: [   2    0   -9    1  ]
# Inequality: [   6   -3    1    3  ]
# Rays 6 5
# Vertex: [ 147  308   33  ]/3
# Vertex: [ 196  407   44  ]/4
# Vertex: [  16   37    4  ]/4
# Vertex: [  98  203   20  ]/2
# Vertex: [  12   28    3  ]/3
# Vertex: [ 147  307   30  ]/3

# UNION LBL: Dimension 3

# LATTICE:
# 4 4
#    2    0    0    0 
#    0    0    0    0 
#   52   66    0    0 
#    0    0    0    1 
# POLYHEDRON Dimension:3
#            Constraints:7  Equations:0  Rays:10  Lines:0
# Constraints 7 5
# Inequality: [ 200  264   -1    0  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [ 674  891    0    3  ]
# Inequality: [   1    0    0   -1  ]
# Inequality: [  -1    0    0   36  ]
# Inequality: [-150 -198    1    0  ]
# Inequality: [  -1    0    5    0  ]
# Rays 10 5
# Vertex: [ 990 -749  198  ]/990
# Vertex: [ 891 -674  198  ]/891
# Vertex: [10692 -8089 2178  ]/297
# Vertex: [ 990 -749  198  ]/33
# Vertex: [10692 -8089 2376  ]/297
# Vertex: [3564 -2696  792  ]/99
# Vertex: [1188 -899  264  ]/99
# Vertex: [1320 -999  264  ]/143
# Vertex: [ 440 -333   88  ]/440
# Vertex: [1188 -899  264  ]/1188
# UNION POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:6  Lines:0
# Constraints 5 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [-150 -198    1    0  ]
# Inequality: [ 674  891    0    3  ]
# Inequality: [  -1    0    0   39  ]
# Inequality: [   1    0    0  -38  ]
# Rays 6 5
# Vertex: [33858 -25615 6930  ]/891
# Vertex: [3861 -2921  792  ]/99
# Vertex: [3861 -2921  858  ]/99
# Vertex: [11583 -8762 2574  ]/297
# Vertex: [33858 -25612 7524  ]/891
# Vertex: [33858 -25615 7524  ]/891
# UNION POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:6  Lines:0
# Constraints 5 5
# Inequality: [ 200  264   -1    0  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [-674 -891    0   -4  ]
# Inequality: [  -1    0    0   36  ]
# Inequality: [  -1    0    5    0  ]
# Rays 6 5
# Vertex: [32076 -24268 7128  ]/891
# Vertex: [160380 -121340 32076  ]/4455
# Vertex: [3960 -2997  792  ]/110
# Vertex: [1760 -1332  352  ]/143
# Vertex: [4752 -3596 1056  ]/297
# Vertex: [1188 -899  264  ]/33
# UNION POLYHEDRON Dimension:3
#            Constraints:6  Equations:0  Rays:8  Lines:0
# Constraints 6 5
# Inequality: [ 200  264   -1    0  ]
# Inequality: [   2    0   -9    0  ]
# Inequality: [  -1    0    0   39  ]
# Inequality: [-674 -891    0   -4  ]
# Inequality: [   1    0    0  -38  ]
# Inequality: [  -1    0    5    0  ]
# Rays 8 5
# Vertex: [33858 -25616 7524  ]/891
# Vertex: [169290 -128080 33858  ]/4455
# Vertex: [15795 -11950 3159  ]/405
# Vertex: [3159 -2390  702  ]/81
# Vertex: [15444 -11687 3432  ]/396
# Vertex: [22572 -17081 5016  ]/594
# Vertex: [8360 -6327 1672  ]/220
# Vertex: [17160 -12987 3432  ]/440
# UNION POLYHEDRON Dimension:3
#            Constraints:5  Equations:0  Rays:6  Lines:0
# Constraints 5 5
# Inequality: [   2    0   -9    0  ]
# Inequality: [-150 -198    1    0  ]
# Inequality: [  -1    0    0   50  ]
# Inequality: [ 674  891    0    3  ]
# Inequality: [   1    0    0  -40  ]
# Rays 6 5
# Vertex: [35640 -26963 7326  ]/891
# Vertex: [44550 -33703 9306  ]/891
# Vertex: [44550 -33703 9900  ]/891
# Vertex: [44550 -33700 9900  ]/891
# Vertex: [35640 -26960 7920  ]/891
# Vertex: [35640 -26963 7920  ]/891


# print("len(YP) =", len(YP))
# if len(YP - diff2) == 0:
#   print("YP is in diff2")
