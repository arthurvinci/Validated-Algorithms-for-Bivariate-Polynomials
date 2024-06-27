################################################################################
#
#   Objects definition
#
################################################################################

using AbstractAlgebra
using Nemo
using ValidatedHensel

# Q[[X]][Y]/(x^10) for exact computation
QX, x = power_series_ring(AbstractAlgebra.QQ, 10, "x"; model=:capped_absolute)
QXY, y = polynomial_ring(QX,"y")

# R[[X]][Y]/(x^10) with 64-bit of precision for approximate computation
RX, _ = power_series_ring(RDF, 10, "x"; model=:capped_absolute)
RXY, _ = polynomial_ring(RX,"y")

################################################################################
#
#   Example: F = (y^2 -(25//10*x^3 -5*x^2 + 845//100*x -576//105))*(y^3 - (x^3 -114//206* x^2 + 30945//109364*x -1))
#            G = y^3 + (9 + 5//4*x + 5*x^2 + 5//3*x^3 + 3//10*x^4)*y^2 + 2*y + 1
#
################################################################################

############################## Exact computation ###############################
Fe = (y^2 -(25//10*x^3 -5*x^2 + 845//100*x -576//105))*(y^3 - (x^3 -114//206* x^2 + 30945//109364*x -1))
Ge = y^3 + (9 + 5//4*x + 5*x^2 + 5//3*x^3 + 3//10*x^4)*y^2 + 2*y + 1

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
l=10

#computation
@time (Qf, Rf, r_Q, r_R) = val_fast_div_rem(Ff, Gf, l, rho, tau)

# Sanity check
Ff - Gf*Qf - Rf


################################# Comparison ###################################
RRDF = typeof(ArbField(128)(1))
RRX, _ = power_series_ring(ArbField(128), 10, "x"; model=:capped_absolute)
RRY, _ = polynomial_ring(RRX, "y")


Qef = to_other_poly(to_other_poly(Qe, RDF, RX, RXY), RRDF, RRX, RRY)
Qff = to_other_poly(Qf, RRDF, RRX, RRY)
actual_Q_bound = mag(biv_norm(Qef - Qff, rho, tau))
given_Q_bound = mag(r_Q)

Ref = to_other_poly(to_other_poly(Re, RDF, RX, RXY), RRDF, RRX, RRY)
Rff = to_other_poly(Rf, RRDF, RRX, RRY)
actual_R_bound = mag(biv_norm(Ref - Rff, rho, tau))
given_R_bound = mag(r_R)