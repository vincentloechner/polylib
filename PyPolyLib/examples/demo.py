from pypolylib import LBL

a = LBL("{(2i, 2i+j) | 1 <= i <= 5, 1 <= j <= i}")
b = LBL("{(i,j) | 0 <= i <= 6, 0 <= j <= 8}")

print(f"a = {a}")
print(f"b = {b}")

print(f"a + b = {a + b}")
print(f"a - b = {a - b}")
print(f"a * b = {a * b}")

print("scan:", end=" ")
for z in a:
  print(z, end= " ")
print()

print(f"set(a) = {set(a)}")

print("a.plot()")
a.plot()

i = a * b
print("plotting (a-i) and i in window1 and (b-i) in window2")
# plot (the border of a) and (a as the two subsets: (a-i) + i) in window1:
a.plot(subplot=True, show_points=False, color="black")
(a-i).plot(subplot=True, color="blue")
i.plot(interactive_update=True, color="cyan")
                                    # leave this window open in background

# plot (the border of b) and (b-i) in window2:
b.plot(subplot=True, show_points=False, color="black")
(b-i).plot(color="green")
