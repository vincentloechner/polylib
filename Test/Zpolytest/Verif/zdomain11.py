import math


def ZP2D(lattice, i_inf, i_sup, j_inf, j_sup, k_inf, k_sup, extra="True"):
  L = set()
  for i in range(math.ceil(eval(i_inf)), math.floor(eval(i_sup))+1):
    for j in range(math.ceil(eval(j_inf)), math.floor(eval(j_sup))+1):
      for k in range(math.ceil(eval(k_inf)), math.floor(eval(k_sup))+1):
        if eval(extra):
          L.add(eval(lattice))
  return L

A = ZP2D("(2*i, 0, 52*i+66*k)", "0", "1000", "i/5", "2*i/9", "(-200*i+j)/264", "(-150*i+j)/198")

print("A contains", len(A), "points")


ZA = set()
# ZA = ZA.union(ZP2D("(2*i, 0, 52*i+66*j)", "0", "1000", "-333*i/440", "(-674*i-4)/891", "0", "0"))
# ZA = ZA.union({(74, 0, 76)}) # <- duplicate!
# ZA = ZA.union({(542, 0, 562)})
# ZA = ZA.union({(0, 0, 0)})
# ZA = ZA.union({(1084, 0, 1124)})
# ZA = ZA.union({(1782, 0, 1848)})
# ZA = ZA.union({(1626, 0, 1686)})
ZA = ZA.union(ZP2D("(2*i, 0, 52*i+66*j)", "0", "1000", "-333*i/440", "(-674*i)/891", "0", "0"))

print("len(ZA) =", len(ZA))

print("A - ZA =", A - ZA)
print("ZA - A =", ZA - A)

