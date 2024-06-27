#######################################################
using AbstractAlgebra
using Nemo
using ValidatedHensel

l = 10
QX, x = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
QXY, y = polynomial_ring(QX,"y")


g = y^3 + (4*x^4 + 5*x)*y^2 + 1
h = y^4 + (3x^2 + 234*x^4+5)*y + x^5

g = y^3 + 5*x*y^2 + 1
h = y^4 + x
F = g*h
g = biv_truncate(g, 2)
h = biv_truncate(h, 2)
s,t = cofactors(g, h, AbstractAlgebra.QQ, QX, QXY)

e = biv_truncate(F -g*h, 2)
g = biv_truncate(g + t*e, 2)
h = biv_truncate(h +s*e, 2)

g = biv_truncate(g, 1, 1)
h = biv_truncate(h, 1, 1)
s = biv_truncate(s, 1, 1)
t = biv_truncate(t, 1, 1)
s*g + t*h


i = 1
j = 1

i *= 2
j *= 2



e = biv_truncate(F -g*h, i, j)
g = biv_truncate(g + t*e, i, j)
h = biv_truncate(h +s*e, i, j)

b = biv_truncate(s*g + t*h -1, i, j)
t = biv_truncate(t - t*b, i, j)
s = biv_truncate(s - s*b, i, j)
