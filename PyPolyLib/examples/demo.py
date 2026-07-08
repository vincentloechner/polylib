from pypolylib import LBL
a = LBL("{(2i, 2i+j) | 1 <= i <= 5, 1 <= j <= i}")
b = LBL("{(i,j)|j<5}")

print(f"a = {a}")
print(f"b = {b}")

print(f"a + b = {a + b}")
print(f"a - b = {a - b}")
print(f"a * b = {a * b}")

print("scan:", end=" ")
for z in a:
  print(z, end= " ")

print(f"set(a) = {set(a)}")

a.plot()

c = LBL("{(i,j) | 0 <= i <= 5, 0 <= j <= 5}")
i = a * c
a1 = a - i
c1 = c - i
(a1+i).plot(interactive_update=True) # let this window close in background
c1.plot()
