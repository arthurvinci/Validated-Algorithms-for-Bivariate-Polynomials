################################################################################
#
#   Objects definition
#
################################################################################

using AbstractAlgebra
using Nemo
using ValidatedHensel

# Q[[X]][Y]/(x^100) for exact computation
QX, x = power_series_ring(AbstractAlgebra.QQ, 10, "x"; model=:capped_absolute)
QXY, y = polynomial_ring(QX,"y")

# R[[X]][Y]/(x^100) with 64-bit of precision for approximate computation
RX, _ = power_series_ring(RDF, 10, "x"; model=:capped_absolute)
RXY, _ = polynomial_ring(RX,"y")

################################################################################
#
#   Example randomly generated 2
#
################################################################################

############################## Exact computation ###############################
Fe = y^100 + (7//2 + O(x^10))*y^63 + (9//8 + 2//3*x^9 + O(x^10))*y^62 + (1//4 + O(x^10))*y^61 + (10//3 + O(x^10))*y^60 + (2 + O(x^10))*y^59 + (2//5 + O(x^10))*y^58 + (2 + O(x^10))*y^57 + (2//3 + O(x^10))*y^56 + (7//5 + O(x^10))*y^55 + y^54 + (2//5 + O(x^10))*y^53 + (8//9 + O(x^10))*y^52 + (5 + O(x^10))*y^51 + (7//2 + O(x^10))*y^49 + (2//3 + O(x^10))*y^48 + y^47 + (6 + O(x^10))*y^46 + (2 + O(x^10))*y^45 + (8//5 + O(x^10))*y^44 + (9//4 + O(x^10))*y^43 + (6//5 + O(x^10))*y^42 + (3//2 + O(x^10))*y^41 + (9//7 + O(x^10))*y^40 + (3//5 + O(x^10))*y^39 + (8//7 + O(x^10))*y^38 + (2//5 + O(x^10))*y^37 + (5//6 + O(x^10))*y^36 + (2//3 + O(x^10))*y^35 + (3//2 + O(x^10))*y^34 + (1//4 + O(x^10))*y^33 + (5//3 + O(x^10))*y^32 + (4//5 + O(x^10))*y^31 + (3 + O(x^10))*y^30 + (7//3 + O(x^10))*y^29 + (6//5 + O(x^10))*y^28 + (1//4 + O(x^10))*y^27 + (3//5 + O(x^10))*y^26 + (8//7 + O(x^10))*y^24 + (3//5 + O(x^10))*y^23 + (4//7 + O(x^10))*y^22 + (7//2 + O(x^10))*y^21 + (5//3 + O(x^10))*y^20 + (7//2 + O(x^10))*y^19 + (7//8 + 5//3*x^3 + 1//2*x^4 + x^5 + 2//9*x^6 + 4//9*x^7 + 9*x^8 + 1//5*x^9 + O(x^10))*y^18 + (3//2 + O(x^10))*y^17 + (2//3 + O(x^10))*y^16 + (7 + 1//3*x^6 + 5//4*x^7 + 2//3*x^8 + 2//3*x^9 + O(x^10))*y^14 + y^12 + (7 + O(x^10))*y^11 + (7//5 + O(x^10))*y^10 + (7 + O(x^10))*y^9 + (4//3 + O(x^10))*y^8 + (1//3 + O(x^10))*y^7 + (1//7 + O(x^10))*y^6 + y^5 + (7//3 + O(x^10))*y^4 + (3//2 + O(x^10))*y^3 + (1//5 + O(x^10))*y^2 + 2 + O(x^10)
Ge = y^75 + (2//3 + O(x^10))*y^61 + (4 + O(x^10))*y^60 + (2 + 3//4*x^3 + 8//7*x^4 + 2//5*x^5 + 7//8*x^6 + 2*x^7 + 7*x^8 + 10//9*x^9 + O(x^10))*y^59 + (5//8 + O(x^10))*y^58 + (6//5 + O(x^10))*y^57 + (9//7 + O(x^10))*y^56 + (9//2 + O(x^10))*y^55 + (3 + O(x^10))*y^54 + (2//9 + 4*x^7 + 4//5*x^8 + 9//8*x^9 + O(x^10))*y^53 + (1//2 + O(x^10))*y^52 + (5//7 + O(x^10))*y^51 + (9//7 + O(x^10))*y^50 + (1//3 + O(x^10))*y^49 + (3 + O(x^10))*y^48 + y^47 + (3 + O(x^10))*y^46 + (1//8 + O(x^10))*y^45 + (7//2 + O(x^10))*y^44 + y^43 + y^42 + (3//4 + O(x^10))*y^41 + (6//7 + O(x^10))*y^40 + (4 + O(x^10))*y^39 + (2 + 2//7*x^7 + 1//8*x^8 + 7*x^9 + O(x^10))*y^37 + (8//7 + O(x^10))*y^36 + (1//8 + O(x^10))*y^35 + (4//3 + O(x^10))*y^34 + y^33 + (10//3 + O(x^10))*y^32 + (5//6 + O(x^10))*y^31 + (2 + O(x^10))*y^30 + (7//3 + O(x^10))*y^29 + (3//5 + O(x^10))*y^28 + (5//3 + O(x^10))*y^27 + (3//8 + O(x^10))*y^26 + (2//3 + O(x^10))*y^25 + (7//8 + O(x^10))*y^24 + (2 + O(x^10))*y^23 + (5//2 + O(x^10))*y^22 + (1//7 + O(x^10))*y^21 + (10//3 + O(x^10))*y^20 + (1//5 + O(x^10))*y^19 + (1//3 + O(x^10))*y^18 + (3//2 + O(x^10))*y^17 + (3//2 + O(x^10))*y^16 + y^15 + (5//9 + O(x^10))*y^14 + (2 + O(x^10))*y^13 + (10//9 + O(x^10))*y^12 + (3//2 + 1//9*x^6 + 4//5*x^7 + 3//2*x^8 + 5//7*x^9 + O(x^10))*y^10 + (7//6 + O(x^10))*y^9 + (2//3 + O(x^10))*y^8 + (1//4 + 10//9*x^5 + 7//9*x^6 + 5*x^9 + O(x^10))*y^7 + (4//3 + O(x^10))*y^6 + (8//9 + O(x^10))*y^4 + (2//3 + O(x^10))*y^3 + (2//3 + O(x^10))*y^2 + (3//5 + 4//7*x^4 + x^5 + 4//3*x^6 + 3//5*x^7 + x^8 + 8//3*x^9 + O(x^10))*y


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