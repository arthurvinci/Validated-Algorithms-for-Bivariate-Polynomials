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
#   Random example 1 
#
################################################################################

########################### Exact computation ############################
ge = y^10 + (4//7 + 5*x^7 + 2//5*x^8 + 3//7*x^9)*y^6 + (5//6 + 10//3*x^3 + 9//8*x^4 + 3//7*x^5 + 8//7*x^7 + 1//9*x^8 + 5//3*x^9)*y^5 + (1//2 + 3//10*x^2 + 6*x^3 + 1//5*x^4 + 3//2*x^5 + 5*x^6 + 3*x^7 + 2//3*x^8)*y^4 + (4//3 + 4*x^4 + 5//2*x^5 + 1//5*x^6 + 8//7*x^7 + 3*x^8 + 7//4*x^9)*y^3 + (3//2)*y^2 + (1//5 + 3//10*x^3 + 3//5*x^4 + 3//2*x^5 + 2*x^6 + 5//3*x^7 + x^8 + 5//9*x^9)*y + 1//2 + 4*x^7 + 3//2*x^8 + 2*x^9
he = y^10 + (7//3 + 7//8*x^7 + 7*x^8 + 8*x^9)*y^8 + (8 + 9//8*x^7 + 2//3*x^8 + 3//4*x^9)*y^7 + (3//4 + 3*x^3 + 5//3*x^4 + 10*x^5 + 3//2*x^7 + 4*x^8 + 3//5*x^9)*y^6 + (4//9 + 2//7*x + 1//10*x^2 + x^3 + 9//8*x^4 + 4//9*x^5 + 7//10*x^6 + 3//5*x^7 + 9//7*x^8 + 2*x^9)*y^5 + (8//5 + 1//2*x^2 + 10//3*x^3 + x^4 + 1//3*x^5 + 2//7*x^6 + 5//7*x^7 + 5//2*x^8 + 5*x^9)*y^4 + (2//5 + 9//8*x^3 + 5//2*x^4 + x^5 + 7//6*x^6 + 5//8*x^7 + 3*x^8 + 7//4*x^9)*y^3 + (10 + 1//2*x^3 + 9//4*x^4 + 3//10*x^5 + 3*x^6 + 7//2*x^7 + 6//7*x^8 + 8*x^9)*y^2 + (3//2 + 6//7*x^2 + 2*x^3 + 6*x^5 + 3//2*x^6 + 2*x^7 + 5//2*x^8 + 8//9*x^9)*y + 5//2 + 4//5*x^6 + 4//7*x^7 + 5//4*x^8 + 9//10*x^9
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

# Sanity checks
sf*gf + tf*hf
biv_truncate(Ff - gf*hf, 1)

# Choice of converging radiuses
rho = 1.0
tau = 1.0
l=10

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