################################################################################
#
#   Objects definition
#
################################################################################

using AbstractAlgebra
using Nemo
using ValidatedHensel

# Q[[X]][Y]/(x^l) for exact computation
l=10
QX, x = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
QXY, y = polynomial_ring(QX,"y")

# R[[X]][Y]/(x^l) with 64-bit of precision for approximate computation
RX, _ = power_series_ring(RDF, l, "x"; model=:capped_absolute)
RXY, _ = polynomial_ring(RX,"y")

################################################################################
#
#   Example: F = y^3+(x^2+4*x)^3 + 2/3 
#            G = y + x^4 + 3x^7 + 35/8
#
################################################################################

############################## Exact computation ###############################
Fe = y^3+(x^2+4*x)^3 + 2//3 
Ge = y + x^4 + 3x^7 + 35//8


# computation
@time (Qe,Re) = fast_div_rem(Fe, Ge, 10)

# Sanity check
Fe - Ge*Qe - Re

############################## Approximate computation ###############################
Ff = to_other_poly(Fe, RDF, RX, RXY)
Gf = to_other_poly(Ge, RDF, RX, RXY)

# Choice of converging radiuses
rho = 1.0
tau = 1.0

#computation
@time (Qf, Rf, r_Q, r_R) = val_fast_div_rem(Ff, Gf, l, rho, tau)
r_Q
r_R
# Sanity check
Ff - Gf*Qf - Rf


################################# Comparison ###################################
Qef = to_other_poly(Qe, RDF, RX, RXY)
actual_Q_bound = biv_norm(Qef - Qf, rho, tau)
given_Q_bound = mag(r_Q)

Ref = to_other_poly(Re, RDF, RX, RXY)
actual_R_bound = biv_norm(Ref - Rf, rho, tau)
given_R_bound = mag(r_R)