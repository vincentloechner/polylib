from pypolylib import LBL

# infinite one:
# a = LBL("{(i, j, k) | 1 <= i <= 5, 1 <= j}")
# a.plot()

a = LBL("{(i, 2j, k) | -2 <= i <= 3, 0 <= j <= 4, 0 <= k <= 2}")
b = LBL("{(i, j, k) |  0 <= i <= 6, 0 <= j <= i, 0 <= k <= i, k<5}")

# diff + the other diff + inter in different colors:
(a-b).plot(subplot=True, color="lightblue")
(b-a).plot(subplot=True, color="yellow")
(a*b).plot(color="darkgreen")
