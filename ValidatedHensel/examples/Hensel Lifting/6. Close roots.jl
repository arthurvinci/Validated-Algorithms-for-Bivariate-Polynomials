using AbstractAlgebra
using Nemo
using ValidatedHensel

################################################################################
#
#   Example where g and h are coprime because of rounding errors
#                 g = (y-eps)^m
#                 h = (y+eps)^m
#   
#                 F = g*h + y^(2m)*x^2
#
################################################################################

############################## Parameters choice ###############################
l = 60
m = 1
eps = 1/10


# Q[[X]][Y]/(x^10) for exact computation
QX, x = power_series_ring(AbstractAlgebra.QQ, 100, "x"; model=:capped_absolute)
QXY, y = polynomial_ring(QX,"y")

# R[[X]][Y]/(x^10) with 64-bit of precision for approximate computation
RX, _ = power_series_ring(RDF, 100, "x"; model=:capped_absolute)
RXY, _ = polynomial_ring(RX,"y")

############################## Exact computation ###############################
ge = (y-eps)^m
he = (y+eps)^m
Fe = ge*he + x^2*y^m

se,te = cofactors(ge, he, AbstractAlgebra.QQ, QX, QXY)

# Sanity checks
se*ge + te*he
biv_truncate(Fe - ge*he, 1)


# computation
@time (Ge, He, Se, Te) = hensel_lifting(Fe, ge, he, se, te, l)

# Sanity checks
biv_truncate(Fe - Ge*He, l) 
biv_truncate(Se*Ge + Te*He, l)


########################### Approximate computation ############################
g = to_other_poly(ge, RDF, RX, RXY)
h = to_other_poly(he, RDF, RX, RXY)
F = to_other_poly(Fe, RDF, RX, RXY)
s = to_other_poly(se, RDF, RX, RXY)
t = to_other_poly(te, RDF, RX, RXY)


# Sanity checks
sf*gf + tf*hf
biv_truncate(Ff - gf*hf, 1)

# Choice of converging radiuses
rho = 0.2
tau = 1.0

# computation
@time (Gf, Hf, r) = val_hensel_lifting(Ff, gf, hf, sf, tf, l, rho, tau)



# Sanity check
biv_norm(Ff - Gf*Hf, rho, tau) 

################################# Comparison ###################################
Gef = to_other_poly(Ge, RDF, RX, RXY)
actual_G_bound = biv_norm(Gef - Gf, rho, tau)
given_G_bound = mag(r)


Hef = to_other_poly(He, RDF, RX, RXY)
actual_H_bound = biv_norm(Hef - Hf, rho, tau)
given_H_bound = mag(r)