################################################################################
#
#   Objects definition
#
################################################################################

using AbstractAlgebra
using Nemo
using ValidatedHensel

# Q[[X]][Y]/(x^10) for exact computation
QX, xe = power_series_ring(AbstractAlgebra.QQ, 10, "x"; model=:capped_absolute)
QXY, ye = polynomial_ring(QX,"y")

# R[[X]][Y]/(x^10) with 64-bit of precision for approximate computation
RX, xf = power_series_ring(RDF, 10, "x"; model=:capped_absolute)
RXY, yf = polynomial_ring(RX,"y")

################################################################################
#
#   Example: (y -(2.5*x^3 -5.123*x^2 + 8.45*x -5.765))*(y - (x^3 -3.14* x^2 + 3.0945*x -1.154))
#
################################################################################

########################### Exact computation ############################
ge = ye - ((25//10)*xe^3 - (5123//1000)*xe^2 + (845//100)*xe - (5765//1000))
he = ye - (xe^3 - (314//100)*xe^2 + (30945//10000)*xe - (1154//1000))
Fe = ge*he

ge = biv_truncate(ge, 1)
he = biv_truncate(he, 1)
se,te = cofactors(ge, he, AbstractAlgebra.QQ, QX, QXY)

# Sanity checks
se*ge + te*he
biv_truncate(Fe - ge*he, 1)

# computation
@time (Ge, He, Se, Te) = hensel_lifting(Fe, ge, he, se, te, 10)

# Sanity checks
Fe - Ge*He 
Se*Ge + Te*He

########################### Approximate computation ############################
gf = to_other_poly(ge, RDF, RX, RXY)
hf = to_other_poly(he, RDF, RX, RXY)
Ff = to_other_poly(Fe, RDF, RX, RXY)
sf = to_other_poly(se, RDF, RX, RXY)
tf = to_other_poly(te, RDF, RX, RXY)

# Sanity check
sf*gf + tf*hf

# Choice of converging radiuses
rho = 1.0
tau = 1.0

# computation
@time (Gf, Hf, r) = val_hensel_lifting(Ff, gf, hf, sf, tf, 10, rho, tau)

# Sanity check
Ff - Gf*Hf 


################################# Comparison ###################################
Gef = to_other_poly(Ge, RDF, RX, RXY)
actual_G_bound = biv_norm(Gef - Gf, rho, tau)
given_G_bound = mag(r)


Hef = to_other_poly(He, RDF, RX, RXY)
actual_H_bound = biv_norm(Hef - Hf, rho, tau)
given_H_bound = mag(r)