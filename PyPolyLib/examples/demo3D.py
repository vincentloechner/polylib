from pypolylib import LBL

# infinite one:
# a = LBL("{(i, j, k) | 1 <= i <= 5, 1 <= j}")
# a.plot()

# a is a rectangular parallelepiped touching the origin on its j face
# (only even j's)
a = LBL("{(i, 2j, k) | -2 <= i <= 3, 0 <= j <= 4, 0 <= k <= 2}")

# b is a dense cone from the origin
b = LBL("{(i, j, k) |  0 <= i <= 6, 0 <= j <= i, 0 <= k <= i, k<5}")

# diff + the other diff + inter in different colors:
(a-b).plot(subplot=True, color="lightblue")
(b-a).plot(subplot=True, color="yellow")
# blue+yellow = green
(a*b).plot(color="darkgreen")
