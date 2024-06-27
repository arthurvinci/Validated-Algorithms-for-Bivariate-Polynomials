using AbstractAlgebra
using Nemo
using ValidatedHensel

################################################################################
#
#   Zheng example g = (y-1)^(4m)*(y-4)^m 
#                 h = (y-2)^(3m)*(y-3)^(2m)
#                 
#                 F = g*h + y^(3*m)*x^2
#
################################################################################

############################## Parameters choice ###############################
l = 16
m = 1

# Q[[X]][Y]/(x^10) for exact computation
QX, x = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
QXY, y = polynomial_ring(QX,"y")

# R[[X]][Y]/(x^10) with 64-bit of precision for approximate computation
RX, _ = power_series_ring(RDF, l, "x"; model=:capped_absolute)
RXY, _ = polynomial_ring(RX,"y")


############################## Exact computation ###############################
ge = (y-1)^(4m)*(y-4)^m
he = (y-2)^(3m)*(y-3)^(2m)
Fe = ge*he + y^(3*m)*x^2

se,te = cofactors(ge, he, AbstractAlgebra.QQ, QX, QXY)

# Sanity checks
se*ge + te*he
biv_truncate(Fe - ge*he, 1)

# computation
@time (Ge, He, Se, Te) = hensel_lifting(Fe, ge, he, se, te, l)
ge
Ge

# Sanity checks
biv_truncate(Fe - Ge*He, l) 
biv_truncate(Se*Ge + Te*He, l)

########################### Approximate computation ############################
gf = to_other_poly(ge, RDF, RX, RXY)
hf = to_other_poly(he, RDF, RX, RXY)
Ff = to_other_poly(Fe, RDF, RX, RXY)
sf = to_other_poly(se, RDF, RX, RXY)
tf = to_other_poly(te, RDF, RX, RXY)


# Sanity checks
sf*gf + tf*hf
biv_truncate(Ff - gf*hf, 1)

# Choice of converging radiuses
rho = 0.005
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