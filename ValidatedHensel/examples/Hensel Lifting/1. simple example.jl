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

# R[[X]][Y]/(x^l) with 64-bit of precision for approximate computation
RX, _ = power_series_ring(RDF, 10, "x"; model=:capped_absolute)
RXY, _ = polynomial_ring(RX,"y")


################################################################################
#
#   Example: (y^15-1)*(y^15+2)+(x^2+4*x)^3
#
################################################################################


############################## Exact computation ###############################
ge = ye^15 - 1
he = ye^15+ 2
Fe = ge*he+(xe^2+4*xe)^3
se,te = cofactors(ge, he, AbstractAlgebra.QQ, QX, QXY)

# Sanity check
se*ge + te*he

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
l = 10

# computation
(Gf, Hf, r) = val_hensel_lifting(Ff, gf, hf, sf, tf, 10, rho, tau)

# Sanity checks
Ff - Gf*Hf 

################################# Comparison ###################################
Gef = to_other_poly(Ge, RDF, RX, RXY)
actual_G_bound = biv_norm(Gef - Gf, rho, tau)
given_G_bound = mag(r)


Hef = to_other_poly(He, RDF, RX, RXY)
actual_H_bound = biv_norm(Hef - Hf, rho, tau)
given_H_bound = mag(r)