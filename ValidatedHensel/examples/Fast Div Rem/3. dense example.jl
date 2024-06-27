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
#   Example randomly generated 
#
################################################################################

############################## Exact computation ###############################
Fe = y^10 + (7//8 + 1//2*x^7 + 2*x^8 + 7//9*x^9 + O(x^10))*y^9 + (3 + 5//4*x^4 + x^5 + 1//4*x^6 + 5*x^7 + 3//5*x^8 + 4//9*x^9 + O(x^10))*y^8 + (1 + 9//5*x^3 + 4//5*x^4 + 8//7*x^5 + 10//3*x^6 + 1//3*x^7 + 1//3*x^8 + 1//2*x^9 + O(x^10))*y^6 + y^5 + (1//4 + O(x^10))*y^4 + (2//7 + x + 2//3*x^2 + 10//7*x^3 + 8//5*x^5 + 5*x^6 + 7//4*x^7 + 3//7*x^8 + 3*x^9 + O(x^10))*y^3 + (1//5 + 5//9*x^2 + 2//5*x^4 + x^5 + 8//9*x^7 + 7*x^8 + 7//9*x^9 + O(x^10))*y^2 + y + 1 + 4//5*x^9 + O(x^10)
Ge = y^5 + (3//7 + 2//3*x^2 + 1//6*x^3 + 2//3*x^4 + 1//3*x^5 + 1//2*x^6 + 1//6*x^7 + 2//5*x^8 + 4//3*x^9 + O(x^10))*y^4 + (2 + x^7 + 2*x^8 + 2//9*x^9 + O(x^10))*y^3 + (5//6 + 7//10*x^2 + 2*x^3 + 5//3*x^4 + 3//7*x^5 + 1//4*x^6 + 1//2*x^7 + 1//9*x^8 + 3//4*x^9 + O(x^10))*y^2 + (5//2 + O(x^10))*y + 2//5 + O(x^10)

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