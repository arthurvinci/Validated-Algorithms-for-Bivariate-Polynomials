using AbstractAlgebra
using Nemo
using ValidatedHensel

CFX, XF = power_series_ring(CFloat, 10, "x"; model=:capped_absolute)
CFXY, YF = polynomial_ring(CFX, "y")
rho, tau = 0.5, 1

F = YF^4 + 4*XF^3*YF^3 + 209*XF^3*YF^2 + XF^6
G = YF^2 + (45XF^5 + 32XF^4)*YF + 34*XF^9 + 134.567

l = 10
k = 15

K = parent(G)
m = F.length - 1
n = G.length - 1

if m < n
    return (K(0), F, 0, 0)
end

FB = to_ball(F, l)
GB = to_ball(G, l)

(H_tilda, r_h) = val_fast_inv(reverse(G), l, m-n+1, rho, tau)
H_tildaB = to_ball(H_tilda, l)

Q_tildaB = reverse(biv_mullow(H_tildaB, reverse(FB), l, m-n+1))
Q_tilda = from_ball(Q_tildaB, l)
r = mag(ball_norm(Q_tildaB - to_ball(Q_tilda, l), rho, tau))
r_Q = mag(ball_norm(FB, rho, tau)*r_h) + r

R_tildaB = biv_truncate(FB, l, n) - biv_mullow(Q_tildaB, GB, l, n)
R_tilda = from_ball(biv_truncate(FB, l, n) - biv_mullow(Q_tildaB, GB, l, n), l)
r = mag(ball_norm(R_tildaB - to_ball(R_tilda, l), rho, tau))
r_R = mag(ball_norm(GB, rho, tau)*r_Q) + r
return(Q_tilda, R_tilda, r_Q, r_R)

F - G*Q_tilda - R_tilda 


################################################################################
#
#   Tests Hensel Exact
#
################################################################################

RX, x = power_series_ring(RDF, 100, "x"; model=:capped_absolute)
RXY, y = polynomial_ring(RX,"y")
F = (y^15-1)*(y^15+2)+(x^2+4*x)^3
g = y^15-1
h = y^15+2
_, s,t = gcdx(g,h)
s*g + t*h 

l = 10
(G, H, S, T) = hensel_lifting(F, g, h, s, t, l)
biv_truncate(F - G*H, l)
S*G + T*H


F = (3*y - (x-2)^10)*(y -(x^2+4)^10) 
g = 3*y - 2^10
h = y - 4^10
_, s, t = gcdx(g,h)
s*g + t*h

l= 20
(G, H, S, T) = hensel_lifting(F, g, h, s, t, l)
biv_truncate(F - G*H, l, 30)
F
G
H
F.length - 1
G.length - 1
H.length
biv_truncate(F - G*H, 20, 30)
S*G + T*H

################################################################################
#
#   Tests Ancien Hensel Approché
#
################################################################################


#################################################################################
RX, x = power_series_ring(RDF, 10, "x"; model=:capped_absolute)
RXY, y = polynomial_ring(RX,"y")
rho = 1.0
tau = 1.0

g = y -(2*x^3 -5*x^2 + 8*x -5)
h = y - (x^3 -3* x^2 + 3*x -1)
F = y^2 + 6*y +5 -10*x
##################################################################################

RX, x = power_series_ring(RDF, 10, "x"; model=:capped_absolute)
RXY, y = polynomial_ring(RX,"y")
F = (y^15-1)*(y^15+2)+(x^2+4*x)^3
g = y^15-1
h = y^15+2
_, s,t = gcdx(g,h)
s*g + t*h 
rho = 1.0
tau = 1.0

l = 10

K = parent(F)
k = parent(F(0))
n = h.length -1
m = F.length - n -1

(G, H, S, T) = hensel_lifting(F, g, h, s, t, l)

FB = to_ball(F, l)
GB = to_ball(G, l)
HB = to_ball(H, l)
SB = to_ball(S, l)
TB = to_ball(T, l) 
gamm = mag(ball_norm(GB, rho, tau))
eta = mag(ball_norm(HB, rho, tau))
sigma = mag(ball_norm(SB, rho, tau))
tau = mag(ball_norm(TB, rho, tau))


EB = biv_truncate(FB, l) - biv_mullow(GB, HB, l)
(Q, R, r_Q, r_R) = val_fast_div_rem(from_ball(SB*EB, l), H, l, rho, tau)
Y_poly = K(0)
set_coefficient!(Y_poly, m + 2*n - 2, k(1))
(PHI, _, r_PHI, _) = val_fast_div_rem(Y_poly, H, l, rho, tau)
r_PHI
PSIB = SB*GB + TB*HB - to_ball(K(1), l)

phi = mag(ball_norm(to_ball(PHI,l), rho, tau)) + r_PHI
psi = mag(ball_norm(PSIB, rho, tau))

delta_G = mag(ball_norm(biv_truncate(TB*EB + to_ball(Q, l)*GB, l, m + 1), rho, tau)) + gamm*r_Q
delta_H = mag(ball_norm(to_ball(R, l), rho, tau)) + r_R

# Compute first polynomial's roots
a_1 = 2*(tau + gamm*sigma*phi)
b_1 = psi + gamm*psi*phi - 1
c_1 = delta_G
disc_1 = b_1*b_1 - 4*a_1*c_1
if disc_1 < 0
    error("Could not validate bounds (first polynomial has no real roots)")
end

r_1_min = (-b_1 - sqrt(disc_1))/(2*a_1)
r_1_max = (-b_1 + sqrt(disc_1))/(2*a_1)

# Compute second polynomial's roots
a_2 = 2*(sigma + eta*sigma*phi) 
b_2 = psi + eta*psi*phi -1
c_2 = delta_H
disc_2 = b_2*b_2 - 4*a_2*c_2
if disc_2 < 0
    error("Could not validate bounds (second polynomial has no real roots)")
end

r_2_min = (-b_2 - sqrt(disc_2))/(2*a_2)
r_2_max = (-b_2 + sqrt(disc_2))/(2*a_2) 

if r_1_min > r_2_max || r_2_min > r_1_max 
    error("Could not validate bounds (roots segments don't intersect)")
end

r = max(r_1_min, r_2_min)





################################################################################
#
#   Tests Nouvel Hensel Approché 1 (PLANTE)
#
################################################################################

RX, x = power_series_ring(RDF, 10, "x"; model=:capped_absolute)
RXY, y = polynomial_ring(RX,"y")
F = (y^15-1)*(y^15+2)+(x^2+4*x)^3
g = y^15-1
h = y^15+2
_, s,t = gcdx(g,h)
s*g + t*h 
rho = 0.5
tau = 1.0


K = parent(F)
k = parent(F(0))
n = h.length -1
m = F.length - n -1
l = 10

(G, H, S, T) = hensel_lifting(F, g, h, s, t, l)

biv_truncate(F - G*H, 10)
biv_truncate(S*G + T*H, 10)

FB = to_ball(F, l)
GB = to_ball(G, l)
HB = to_ball(H, l)
SB = to_ball(S, l)
TB = to_ball(T, l) 
norm_G = mag(ball_norm(GB, rho, tau))
norm_H = mag(ball_norm(HB, rho, tau))
norm_S = mag(ball_norm(SB, rho, tau))
norm_T = mag(ball_norm(TB, rho, tau))


EB = biv_truncate(FB, l) - biv_mullow(GB, HB, l)
(Q, R, r_Q, r_R) = val_fast_div_rem(from_ball(SB*EB, l), H, l, rho, tau)

biv_truncate(TB*EB + to_ball(Q, l)*GB, l, m + 1)
delta_G = mag(ball_norm(biv_truncate(TB*EB + to_ball(Q, l)*GB, l, m + 1), rho, tau)) + norm_G*r_Q
delta_H = mag(ball_norm(to_ball(R, l), rho, tau)) + r_R
delta = max(delta_G, delta_H)

# Lipschitz ratio bounding
Y_poly = K(0)
set_coefficient!(Y_poly, m + 2*n - 2, k(1))
(PHI, _, r_PHI, _) = val_fast_div_rem(Y_poly, H, l, rho, tau)
r_PHI

PSIB = SB*GB + TB*HB - to_ball(K(1), l)

mu = max(norm_G, norm_H)
nu = max(norm_T, norm_S)
phi = mag(ball_norm(to_ball(PHI,l), rho, tau)) + r_PHI
psi = mag(ball_norm(PSIB, rho, tau))

# Compute the polynomials' roots
a = 2*nu*(1+mu*phi)
b = psi*(1 + mu*phi) - 1
c = delta
disc = b*b - 4*a*c
if disc < 0
    error("Could not validate bounds (polynomial has no real roots)")
end

r_min = (-b - sqrt(disc))/(2*a)
r_min_2 = -2*c/(b - sqrt(disc))
r_max = (-b + sqrt(disc))/(2*a)
r_max_2 = -2*c/(b + sqrt(disc))
r = r_min
if r_min < 0
    if r_max < 0
        error("Could not validate bounds (polynomial has no positive roots)")
    end
    r = r_max
end

################################################################################
#
#   Tests Nouvel Hensel Approché 2 (RÉSULTAT DOUTEUX)
#
################################################################################
RX, x = power_series_ring(RDF, 10, "x"; model=:capped_absolute)
RXY, y = polynomial_ring(RX,"y")
rho = 0.5
tau = 1.0

g = y -(2*x^3 -5*x^2 + 8*x -5)
h = y - (x^3 -3* x^2 + 3*x -1)
F = g*h

_, s, t = gcdx(g,h)
s*g + t*h

K = parent(F)
k = parent(F(0))
n = h.length -1
m = F.length - n -1

(G, H, S, T) = hensel_lifting(F, g, h, s, t, l)

biv_truncate(F - G*H, 10)
biv_truncate(S*G + T*H, 10)

FB = to_ball(F, l)
GB = to_ball(G, l)
HB = to_ball(H, l)
SB = to_ball(S, l)
TB = to_ball(T, l) 
norm_G = mag(ball_norm(GB, rho, tau))
norm_H = mag(ball_norm(HB, rho, tau))
norm_S = mag(ball_norm(SB, rho, tau))
norm_T = mag(ball_norm(TB, rho, tau))


EB = biv_truncate(FB, l) - biv_mullow(GB, HB, l)
(Q, R, r_Q, r_R) = val_fast_div_rem(from_ball(SB*EB, l), H, l, rho, tau)

biv_truncate(TB*EB + to_ball(Q, l)*GB, l, m + 1)
delta_G = mag(ball_norm(biv_truncate(TB*EB + to_ball(Q, l)*GB, l, m + 1), rho, tau)) + norm_G*r_Q
delta_H = mag(ball_norm(to_ball(R, l), rho, tau)) + r_R
delta = max(delta_G, delta_H)

# Lipschitz ratio bounding
Y_poly = K(0)
set_coefficient!(Y_poly, m + 2*n - 2, k(1))
(PHI, _, r_PHI, _) = val_fast_div_rem(Y_poly, H, l, rho, tau)


PSIB = SB*GB + TB*HB - to_ball(K(1), l)

mu = max(norm_G, norm_H)
nu = max(norm_T, norm_S)
phi = mag(ball_norm(to_ball(PHI,l), rho, tau)) + r_PHI
psi = mag(ball_norm(PSIB, rho, tau))

# Compute the polynomials' roots
a = 2*nu*(1+mu*phi)
b = psi*(1 + mu*phi) - 1
c = delta
disc = b*b - 4*a*c
if disc < 0
    error("Could not validate bounds (polynomial has no real roots)")
end

r_min = (-b - sqrt(disc))/(2*a)
r_max = (-b + sqrt(disc))/(2*a)

r = r_min
if r_min < 0
    if r_max < 0
        error("Could not validate bounds (polynomial has no positive roots)")
    end
    r = r_max
end


################################################################################
#
#   Tests Nouvel Hensel Approché Aléatoire (PLANTE)
#
################################################################################
RX, x = power_series_ring(RDF, 10, "x"; model=:capped_absolute)
RXY, y = polynomial_ring(RX,"y")

g = random_polynomial(RXY, RX, 2, 2)
h = random_polynomial(RXY, RX, 2, 2)
F = biv_mullow(g,h, 10)
g = biv_truncate(g, 1)
h = biv_truncate(h, 1)


_, s,t = gcdx(g,h)
s*g + t*h 
rho = 0.1
tau = 1.0

F - g*h

l=10
K = parent(F)
k = parent(F(0))
n = h.length -1
m = F.length - n -1

(G, H, S, T) = hensel_lifting(F, g, h, s, t, l)

biv_truncate(F - G*H, 10)
biv_truncate(S*G + T*H, 10)

FB = to_ball(F, l)
GB = to_ball(G, l)
HB = to_ball(H, l)
SB = to_ball(S, l)
TB = to_ball(T, l) 
norm_G = mag(ball_norm(GB, rho, tau))
norm_H = mag(ball_norm(HB, rho, tau))
norm_S = mag(ball_norm(SB, rho, tau))
norm_T = mag(ball_norm(TB, rho, tau))


EB = biv_truncate(FB, l) - biv_mullow(GB, HB, l)
(Q, R, r_Q, r_R) = val_fast_div_rem(from_ball(SB*EB, l), H, l, rho, tau)

biv_truncate(TB*EB + to_ball(Q, l)*GB, l, m + 1)
delta_G = mag(ball_norm(biv_truncate(TB*EB + to_ball(Q, l)*GB, l, m + 1), rho, tau)) + norm_G*r_Q
delta_H = mag(ball_norm(to_ball(R, l), rho, tau)) + r_R
delta = max(delta_G, delta_H)

# Lipschitz ratio bounding
Y_poly = K(0)
set_coefficient!(Y_poly, m + 2*n - 2, k(1))
(PHI, _, r_PHI, _) = val_fast_div_rem(Y_poly, H, l, rho, tau)
r_PHI

PSIB = SB*GB + TB*HB - to_ball(K(1), l)

mu = max(norm_G, norm_H)
nu = max(norm_T, norm_S)
phi = mag(ball_norm(to_ball(PHI,l), rho, tau)) + r_PHI
psi = mag(ball_norm(PSIB, rho, tau))

# Compute the polynomials' roots
a = 2*nu*(1+mu*phi)
b = psi*(1 + mu*phi) - 1
c = delta
disc = b*b - 4*a*c
if disc < 0
    error("Could not validate bounds (polynomial has no real roots)")
end

r_min = (-b - sqrt(disc))/(2*a)
r_min_2 = -2*c/(b - sqrt(disc))
r_max = (-b + sqrt(disc))/(2*a)
r_max_2 = -2*c/(b + sqrt(disc))

r = r_min
if r_min < 0
    if r_max < 0
        error("Could not validate bounds (polynomial has no positive roots)")
    end
    r = r_max
end



################################################################################
#
#   Tests Ancien Hensel Approché
#
################################################################################

RX, x = power_series_ring(RDF, 10, "x"; model=:capped_absolute)
RXY, y = polynomial_ring(RX,"y")
F = (y^15-1)*(y^15+2)+(x^2+4*x)^3
g = y^15-1
h = y^15+2
_, s,t = gcdx(g,h)
s*g + t*h 
rho = 0.5
tau = 1.0

K = parent(F)
k = parent(F(0))
n = h.length -1
m = F.length - n -1

(G, H, S, T) = hensel_lifting(F, g, h, s, t, l)

FB = to_ball(F, l)
GB = to_ball(G, l)
HB = to_ball(H, l)
SB = to_ball(S, l)
TB = to_ball(T, l) 
gamm = mag(ball_norm(GB, rho, tau))
eta = mag(ball_norm(HB, rho, tau))
sigma = mag(ball_norm(SB, rho, tau))
tau = mag(ball_norm(TB, rho, tau))


EB = biv_truncate(FB, l) - biv_mullow(GB, HB, l)
(Q, R, r_Q, r_R) = val_fast_div_rem(from_ball(SB*EB, l), H, l, rho, tau)
Y_poly = K(0)
set_coefficient!(Y_poly, m + 2*n - 2, k(1))
(PHI, _, r_PHI, _) = val_fast_div_rem(Y_poly, H, l, rho, tau)
r_PHI
PSIB = SB*GB + TB*HB - to_ball(K(1), l)

phi = mag(ball_norm(to_ball(PHI,l), rho, tau)) + r_PHI
psi = mag(ball_norm(PSIB, rho, tau))

delta_G = mag(ball_norm(biv_truncate(TB*EB + to_ball(Q, l)*GB, l, m + 1), rho, tau)) + gamm*r_Q
delta_H = mag(ball_norm(to_ball(R, l), rho, tau)) + r_R

# Compute first polynomial's roots
a_1 = 2*(tau + gamm*sigma*phi)
b_1 = psi + gamm*psi*phi - 1
c_1 = delta_G
disc_1 = b_1*b_1 - 4*a_1*c_1
if disc_1 < 0
    error("Could not validate bounds (first polynomial has no real roots)")
end

r_1_min = (-b_1 - sqrt(disc_1))/(2*a_1)
r_1_max = (-b_1 + sqrt(disc_1))/(2*a_1)

# Compute second polynomial's roots
a_2 = 2*(sigma + eta*sigma*phi) 
b_2 = psi + eta*psi*phi -1
c_2 = delta_H
disc_2 = b_2*b_2 - 4*a_2*c_2
if disc_2 < 0
    error("Could not validate bounds (second polynomial has no real roots)")
end

r_2_min = (-b_2 - sqrt(disc_2))/(2*a_2)
r_2_max = (-b_2 + sqrt(disc_2))/(2*a_2) 

if r_1_min > r_2_max || r_2_min > r_1_max 
    error("Could not validate bounds (roots segments don't intersect)")
end

r = max(r_1_min, r_2_min)


###################################################################################################

k = 30
l= 20
# Q[[X]][Y]/(x^l) for exact computation
QX, x = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
QXY, y = polynomial_ring(QX,"y")
gi = y^30 + (1//10 + O(x^20))*y + 2//5 + O(x^20)
hi = y^30 + (3//7 + O(x^20))*y^23 + (3 + O(x^20))*y^21 + (3//8 + O(x^20))*y^20 + (2 + O(x^20))*y^19 + (8//5 + O(x^20))*y^18 + y^17 + (2//3 + O(x^20))*y^16 + (5//7 + O(x^20))*y^15 + (3//4 + O(x^20))*y^14 + (4//7 + O(x^20))*y^13 + (9//10 + O(x^20))*y^12 + (1//2 + O(x^20))*y^11 + (7//9 + O(x^20))*y^9 + (3//2 + O(x^20))*y^8 + y^7 + (3//4 + O(x^20))*y^6 + (3//2 + O(x^20))*y^5 + (5//4 + O(x^20))*y^4 + (7//6 + O(x^20))*y^3 + (6//5 + O(x^20))*y + 1//10 + O(x^20)
F = g*h

g = biv_truncate(gi, 1)
h = biv_truncate(hi, 1)
s,t = cofactors(g, h, AbstractAlgebra.QQ, QX, QXY)

(G, H, S, T) = hensel_lifting(F, g, h, s, t, l)

gi == G
hi == H
F - G*H
G*S + T*H


#################################################
k = 100
l = 100
rho = 0.1
tau = 1.0

QX, x = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
QXY, y = polynomial_ring(QX,"y")


A = y^100 + (1//3 + 2//9*x^7 + x^8 + 1//2*x^9 + 7//5*x^10 + 9//4*x^11 + 10//7*x^12 + 6//7*x^13 + 3//2*x^14 + 6*x^18 + 4*x^19 + 6*x^20 + 1//2*x^21 + 9//8*x^22 + 2//5*x^23 + x^24 + 3//10*x^25 + x^26 + 7//5*x^27 + 9//10*x^28 + 7//5*x^29 + 4//7*x^30 + 10//7*x^31 + 9//4*x^32 + 9*x^33 + 8*x^34 + 4//5*x^35 + 7//6*x^36 + 2//9*x^37 + 3//2*x^38 + 8//3*x^39 + x^40 + 5*x^41 + 9//2*x^42 + 3//2*x^43 + 5//2*x^44 + 3//4*x^45 + 3//5*x^46 + 4*x^47 + 5//2*x^48 + x^49 + 2//5*x^50 + 7//2*x^51 + 2*x^52 + x^54 + 9//5*x^55 + 9//2*x^56 + 9*x^57 + 8*x^58 + 2*x^59 + 1//4*x^60 + 7*x^61 + 5//3*x^62 + 7//10*x^63 + 2//5*x^64 + 1//2*x^65 + 3*x^66 + 5//3*x^67 + 1//8*x^68 + x^69 + x^70 + 2*x^72 + 1//2*x^73 + 1//3*x^74 + 10//3*x^75 + 3//2*x^76 + 4*x^77 + 2//7*x^78 + x^79 + 3//5*x^80 + 2//3*x^81 + 1//5*x^82 + 4//3*x^83 + 7//6*x^85 + 5//2*x^86 + 5//2*x^87 + 3//2*x^88 + 10//3*x^89 + 2//5*x^90 + 4//9*x^91 + 3//4*x^92 + 3//4*x^93 + x^94 + 10//3*x^95 + 4//5*x^96 + 2//9*x^97 + 1//5*x^98 + 8//7*x^99 + O(x^100))*y^25 + (3 + 1//2*x^72 + 7*x^73 + x^74 + 1//4*x^75 + 5//4*x^77 + 9//5*x^78 + 5*x^79 + 2//3*x^80 + 5//3*x^81 + 5//4*x^82 + x^83 + 8*x^84 + 6*x^85 + 2*x^86 + 2//5*x^87 + 10*x^88 + 5//9*x^89 + 3//4*x^90 + 1//4*x^91 + 4//3*x^92 + 5//4*x^93 + x^94 + 5//6*x^95 + 4//5*x^96 + x^97 + 1//2*x^98 + 2*x^99 + O(x^100))*y^24 + (1//10 + 5//9*x^94 + x^95 + 3//2*x^96 + 1//7*x^97 + 5//2*x^98 + 6//7*x^99 + O(x^100))*y^23 + (4 + 2//3*x^88 + 8//7*x^89 + x^90 + 8//7*x^91 + 9//4*x^93 + 7//8*x^94 + 3//4*x^95 + 3//2*x^96 + 3//8*x^97 + 3//8*x^98 + 2*x^99 + O(x^100))*y^21 + (1//9 + x^89 + 1//4*x^90 + 7//8*x^92 + x^93 + 5//4*x^95 + 1//2*x^96 + x^98 + 6//5*x^99 + O(x^100))*y^20 + (1 + 2//9*x^92 + 2*x^93 + 7//2*x^94 + 2//5*x^96 + x^97 + 5//4*x^99 + O(x^100))*y^19 + (5//8 + 3*x^38 + 4//3*x^40 + 5*x^41 + 3//5*x^42 + 7*x^43 + 6//5*x^44 + 2//5*x^45 + 5//8*x^46 + 1//4*x^47 + 1//3*x^48 + 1//7*x^49 + 7//4*x^51 + 5//4*x^52 + 1//2*x^53 + 5*x^54 + 7//8*x^55 + 1//9*x^56 + 9*x^57 + 1//2*x^58 + 5//8*x^59 + 1//2*x^60 + 6*x^61 + 2//5*x^62 + 3//5*x^63 + x^64 + x^65 + 7//5*x^66 + 2*x^67 + x^68 + 3//5*x^69 + 7//8*x^70 + 3*x^71 + 5//4*x^72 + 1//3*x^73 + x^74 + 5//2*x^75 + 5//2*x^76 + 3//5*x^77 + 9//5*x^78 + 3//7*x^79 + 2//7*x^80 + 2//9*x^81 + 2*x^82 + 9//7*x^83 + 1//6*x^84 + 9//5*x^85 + x^86 + 2//5*x^87 + x^89 + 2*x^90 + 7*x^91 + 9//7*x^92 + x^93 + 1//3*x^94 + 9*x^95 + 1//2*x^96 + 5//4*x^97 + 5//4*x^98 + 9//10*x^99 + O(x^100))*y^18 + (1//10 + 1//3*x^82 + 7//4*x^83 + 7//9*x^84 + 5//2*x^85 + 2//7*x^86 + 3*x^87 + 9//4*x^88 + 5//8*x^89 + 7//8*x^90 + 5//8*x^91 + 2//5*x^92 + 1//2*x^93 + 4//5*x^94 + 2*x^95 + 3//10*x^96 + 2*x^97 + 1//2*x^98 + 1//2*x^99 + O(x^100))*y^17 + (4//5 + 2//3*x^4 + 8*x^5 + 7//3*x^7 + 9//4*x^8 + 5*x^9 + 2//9*x^10 + 7//9*x^11 + 3//8*x^12 + 2//3*x^14 + 9//7*x^15 + 1//10*x^16 + 5//3*x^17 + 2*x^18 + 3//5*x^19 + 5*x^20 + 8//3*x^21 + 5//3*x^22 + 7//5*x^23 + 2//7*x^24 + 2//5*x^25 + 5//9*x^26 + 5//3*x^27 + 5*x^28 + 1//2*x^29 + 1//10*x^30 + 8*x^31 + 5*x^32 + 7//3*x^33 + 4//5*x^34 + 5*x^35 + 2//7*x^36 + 1//8*x^37 + 9//2*x^38 + 5*x^39 + x^40 + 1//7*x^41 + 2*x^42 + 7//9*x^43 + x^44 + 8//9*x^46 + x^47 + 1//10*x^48 + 3//5*x^49 + 2*x^50 + 9//10*x^51 + 1//4*x^53 + x^54 + 8*x^55 + 8*x^56 + 5//2*x^57 + 3//7*x^58 + 5//8*x^59 + 2*x^60 + 7//5*x^61 + 5//4*x^62 + 2//3*x^63 + 3*x^64 + 1//5*x^65 + 4//5*x^66 + 2//9*x^67 + 5//2*x^68 + 7*x^70 + 1//2*x^71 + 2//3*x^73 + 8//5*x^74 + 1//7*x^75 + 3*x^76 + 1//5*x^77 + 2//9*x^78 + 2//5*x^79 + 1//8*x^80 + 1//10*x^81 + 1//2*x^82 + 1//2*x^83 + 1//5*x^84 + 8//9*x^85 + 7//6*x^86 + 6//7*x^87 + 7//5*x^88 + 5//2*x^89 + 4*x^90 + 7//4*x^91 + 2*x^92 + 1//2*x^93 + 2//9*x^94 + x^95 + 1//2*x^98 + 4//5*x^99 + O(x^100))*y^16 + (6//5 + 4//5*x^15 + 5//4*x^16 + 10//7*x^17 + 2//5*x^18 + 9//4*x^19 + 1//2*x^20 + 1//8*x^22 + 7*x^23 + 7//10*x^26 + 8//9*x^27 + 9//8*x^28 + 5//4*x^29 + 9//7*x^30 + 9//4*x^31 + x^32 + 7//5*x^33 + 9//5*x^34 + 1//2*x^35 + 5*x^36 + 7//5*x^37 + 9*x^38 + 2*x^39 + 1//2*x^40 + 9//10*x^41 + 7//2*x^42 + x^43 + 3*x^44 + 1//6*x^45 + 2*x^47 + 1//3*x^48 + 4//3*x^49 + 2*x^50 + x^52 + 6*x^53 + 6*x^55 + 1//3*x^56 + 3//10*x^57 + x^58 + 8//7*x^59 + 1//5*x^60 + 1//3*x^63 + 4//7*x^64 + 7//9*x^65 + 5//2*x^66 + 5//3*x^67 + 7//5*x^68 + 3//2*x^69 + x^70 + 5//2*x^72 + 2*x^74 + 1//10*x^75 + 8//3*x^76 + 2//3*x^77 + 4//3*x^78 + 5//8*x^79 + x^80 + 5//3*x^81 + 9//5*x^82 + 2//9*x^83 + 1//10*x^84 + 2//3*x^85 + 9//10*x^86 + 5//3*x^87 + 9//5*x^88 + 1//3*x^89 + x^90 + 9//5*x^91 + 1//7*x^92 + 1//5*x^93 + 9//2*x^94 + 7//8*x^95 + 8//9*x^96 + 2*x^97 + 2*x^98 + 8//3*x^99 + O(x^100))*y^15 + (10//9 + 5//3*x^80 + 9//8*x^81 + 7*x^82 + 2*x^83 + 3*x^84 + x^85 + 1//10*x^86 + 3*x^87 + 3*x^88 + 5//3*x^91 + 3//10*x^93 + 4//5*x^94 + 5//8*x^95 + 3*x^96 + 9*x^97 + 5//2*x^98 + 5//2*x^99 + O(x^100))*y^14 + (1 + 1//3*x^43 + 6//7*x^44 + x^45 + 5//3*x^46 + x^47 + 1//10*x^48 + 1//3*x^49 + 3//2*x^51 + 9*x^52 + 1//9*x^53 + 5//6*x^54 + 7//10*x^55 + 5//2*x^56 + x^57 + 3*x^58 + 5//8*x^59 + 5//4*x^60 + 7*x^61 + 6//7*x^62 + 7//4*x^63 + 7//6*x^64 + x^65 + 10*x^66 + 9//10*x^67 + 9//7*x^68 + 7//9*x^69 + 6//7*x^70 + 5//2*x^71 + 3//2*x^72 + 3*x^73 + x^74 + x^75 + 9//5*x^76 + 5//8*x^77 + 6//7*x^78 + 4//9*x^79 + 2//3*x^80 + 5//4*x^81 + 1//4*x^82 + 5//3*x^83 + 1//7*x^84 + 10*x^85 + x^86 + 3*x^87 + 3//8*x^88 + x^89 + 5*x^90 + 7//6*x^92 + 8//5*x^93 + 6//5*x^94 + 3*x^95 + 1//2*x^96 + 3//2*x^97 + 5//3*x^98 + 2//9*x^99 + O(x^100))*y^13 + (2 + 1//5*x^72 + 7//5*x^73 + 1//2*x^74 + 4//3*x^75 + 2//7*x^76 + 9//2*x^77 + 5//9*x^78 + 7//3*x^79 + 10//3*x^80 + 7//10*x^81 + 3//4*x^82 + 5//2*x^83 + 1//10*x^84 + x^85 + 6//7*x^86 + 3*x^87 + 7//9*x^88 + 7//6*x^89 + 1//5*x^90 + 1//4*x^91 + 5//9*x^92 + 2*x^93 + x^94 + x^95 + x^96 + 3*x^97 + 5//3*x^98 + 9//4*x^99 + O(x^100))*y^12 + (4//5 + 10*x^13 + 2*x^14 + 1//10*x^15 + 3//4*x^16 + 3*x^17 + x^19 + 10//9*x^20 + 1//10*x^21 + 2*x^22 + 6//7*x^23 + 2*x^25 + 5*x^26 + 1//5*x^27 + 3//4*x^28 + 4//5*x^29 + 1//2*x^30 + 1//9*x^31 + x^33 + 1//3*x^35 + x^37 + 5//9*x^39 + 8*x^40 + 3//7*x^41 + 8*x^43 + 1//6*x^44 + 1//7*x^45 + 5//9*x^46 + 3//10*x^47 + x^48 + 1//3*x^50 + 5*x^51 + 3*x^52 + 10//9*x^53 + 9//5*x^55 + 8//3*x^56 + 3//7*x^57 + 1//8*x^58 + 5//6*x^59 + 5//3*x^61 + 6*x^62 + 2//3*x^63 + 1//9*x^64 + 5//6*x^66 + 7*x^67 + 4//7*x^68 + 9//4*x^70 + 7//3*x^71 + x^72 + 3//2*x^73 + 1//9*x^74 + 6//7*x^75 + 5//8*x^76 + 9//2*x^77 + 6*x^79 + 2//3*x^80 + 5//4*x^81 + x^82 + 2*x^84 + 2//3*x^85 + 1//9*x^86 + 7//8*x^87 + 10*x^88 + 2//3*x^89 + 2//3*x^90 + x^91 + 2//3*x^92 + 1//3*x^93 + 3//7*x^94 + 1//8*x^95 + 4*x^96 + 3//2*x^97 + 8//3*x^98 + 4//5*x^99 + O(x^100))*y^11 + (7//3 + 3//4*x^2 + x^4 + 3//7*x^5 + 3//7*x^6 + 1//3*x^7 + 1//2*x^8 + 1//2*x^9 + 2*x^10 + 3//2*x^11 + 2*x^12 + 9//7*x^13 + 5//8*x^14 + 5//4*x^15 + 4//3*x^16 + 7//5*x^17 + 2//5*x^18 + 4//3*x^19 + 9//7*x^20 + 9*x^21 + 2*x^22 + 9*x^23 + 3//5*x^24 + 7//2*x^25 + x^26 + 10*x^29 + 10//7*x^30 + 7//2*x^31 + 1//2*x^32 + 9*x^33 + 7//8*x^34 + 3//10*x^35 + 1//7*x^36 + 1//2*x^37 + 10//3*x^38 + x^39 + 3//2*x^40 + 3//4*x^41 + x^42 + 2*x^43 + x^44 + 4//3*x^45 + 3//5*x^46 + x^47 + 7*x^48 + 9//2*x^49 + 7//2*x^50 + x^51 + 1//10*x^52 + 5//2*x^53 + 1//5*x^54 + 1//3*x^55 + 1//3*x^56 + x^57 + 2//5*x^58 + 4//9*x^59 + 1//5*x^60 + 2*x^61 + 1//9*x^62 + 7//8*x^63 + x^64 + 1//4*x^65 + 3//4*x^66 + 1//6*x^67 + 2//7*x^68 + 7//9*x^69 + 3//2*x^70 + 3//2*x^71 + 4//3*x^72 + 1//9*x^73 + x^75 + 1//2*x^76 + 10*x^77 + 3*x^78 + 1//5*x^79 + 10//9*x^80 + 2*x^81 + x^82 + 7//10*x^84 + 6*x^85 + 9//2*x^86 + 3//4*x^87 + 4//3*x^88 + 1//6*x^89 + 3//4*x^90 + 3//4*x^91 + 10*x^92 + 4//5*x^93 + x^94 + 5//6*x^96 + 2//3*x^97 + 4//9*x^98 + O(x^100))*y^10 + (1 + 7//10*x^77 + 2//7*x^78 + 3//5*x^79 + 2*x^80 + 1//3*x^81 + x^82 + 3//4*x^84 + 9//8*x^85 + 2//3*x^86 + 1//2*x^87 + 3//8*x^88 + 1//5*x^89 + 2//3*x^90 + 6//7*x^91 + 8//7*x^92 + x^94 + x^95 + 5*x^96 + 7//9*x^97 + 4//9*x^98 + x^99 + O(x^100))*y^9 + (6//7 + 2//7*x^35 + 7//6*x^36 + 7//6*x^37 + 5//9*x^39 + 9//4*x^40 + 7//5*x^42 + 1//9*x^43 + 5//3*x^44 + 3//10*x^45 + 9//8*x^46 + 1//6*x^47 + 2*x^48 + 2*x^50 + 1//3*x^51 + 1//4*x^52 + 6*x^53 + 5//4*x^54 + 1//4*x^55 + 5//4*x^56 + 1//2*x^57 + 9//2*x^58 + 5//7*x^59 + 3//7*x^60 + 5//4*x^61 + 1//5*x^62 + 2//3*x^63 + 9//4*x^64 + 9//7*x^66 + 6*x^67 + x^68 + 2*x^69 + 5//8*x^70 + 8//5*x^71 + 5*x^72 + 2//5*x^73 + x^74 + 3//8*x^75 + 4//9*x^76 + 3//2*x^77 + 8//7*x^78 + 3//10*x^79 + 5//2*x^80 + 4//5*x^81 + 5//6*x^82 + 10//3*x^83 + 1//2*x^84 + 1//4*x^85 + 2*x^86 + 3//2*x^87 + 7//2*x^88 + 5//4*x^89 + 1//5*x^90 + 2//3*x^91 + 1//2*x^92 + 5//4*x^93 + 5//8*x^94 + x^95 + 8//7*x^96 + 4//3*x^97 + 3//5*x^98 + 1//2*x^99 + O(x^100))*y^8 + (5//7 + 2*x^29 + 2//9*x^30 + 3//5*x^31 + 2//5*x^32 + 5//2*x^34 + 10*x^35 + 10*x^36 + 10//3*x^37 + 4//9*x^38 + 2//3*x^41 + 5*x^43 + 1//2*x^44 + 7//3*x^45 + 8//3*x^46 + 10*x^47 + 3//4*x^48 + 5//7*x^49 + 2*x^50 + x^51 + 1//6*x^53 + 10//3*x^54 + 8//5*x^55 + 1//2*x^56 + 10*x^57 + x^58 + 1//7*x^59 + 3//5*x^60 + 3*x^61 + 3//5*x^62 + 6*x^63 + 9//8*x^64 + 8//7*x^65 + 9//4*x^66 + 1//8*x^67 + 5//4*x^68 + 1//3*x^69 + 1//2*x^70 + 3//8*x^71 + 1//2*x^72 + 2//7*x^73 + 7//3*x^74 + x^75 + 5//2*x^76 + 2*x^77 + 1//2*x^78 + 5//8*x^79 + 10//7*x^81 + 2*x^82 + 2*x^83 + 9*x^84 + 4//5*x^85 + 4//3*x^86 + 5//8*x^87 + x^88 + 1//10*x^89 + x^90 + 5//8*x^91 + 3//5*x^92 + 1//9*x^93 + 9//2*x^94 + 4//9*x^95 + 9//2*x^97 + 5*x^98 + O(x^100))*y^7 + (1 + x^96 + 3//10*x^97 + 1//3*x^98 + 4*x^99 + O(x^100))*y^6 + (2//3 + 4//5*x^59 + 1//4*x^60 + 1//2*x^61 + 3//2*x^62 + 6//7*x^63 + 2//5*x^65 + 7//5*x^66 + 2//5*x^67 + 3//2*x^68 + 1//3*x^69 + 3//7*x^70 + 5//2*x^71 + 7//5*x^72 + 5//2*x^73 + 3//2*x^74 + 9//8*x^76 + 2*x^77 + 4//5*x^78 + 2*x^81 + 1//2*x^82 + 8//5*x^83 + 7//4*x^84 + 8//5*x^85 + 9*x^86 + 4*x^87 + 2//3*x^88 + x^89 + 3//2*x^90 + 1//5*x^91 + 7//5*x^92 + 2//3*x^94 + 2//5*x^95 + 1//9*x^96 + 1//2*x^97 + 2*x^98 + x^99 + O(x^100))*y^5 + (8//3 + x^81 + 3//8*x^82 + 7*x^83 + 8//9*x^85 + 8//9*x^86 + 5//2*x^87 + 1//4*x^90 + 3//2*x^91 + 5*x^92 + x^93 + 3//5*x^94 + x^95 + x^96 + x^97 + 5//2*x^98 + x^99 + O(x^100))*y^4 + (2 + x^87 + 2//5*x^88 + 5//2*x^90 + 7*x^92 + 9//2*x^94 + 4*x^95 + 7//3*x^96 + 2//5*x^97 + 4//3*x^98 + x^99 + O(x^100))*y^3 + (1//3 + 8//5*x^52 + 1//2*x^53 + x^56 + 4//7*x^57 + 8//5*x^58 + 2//7*x^59 + 5*x^60 + x^62 + 5//4*x^63 + 1//3*x^64 + 7//8*x^65 + x^66 + 9//2*x^67 + 8//9*x^69 + 6//5*x^70 + 7//2*x^71 + 5*x^72 + 9//5*x^73 + x^74 + 2//3*x^75 + x^76 + 2*x^77 + 3//7*x^78 + 8//3*x^79 + x^80 + 2*x^81 + 7//10*x^82 + 3*x^83 + 2//3*x^84 + 5//8*x^85 + 8//3*x^86 + 8//3*x^87 + 1//6*x^88 + 2//7*x^89 + 7//9*x^91 + 2//5*x^92 + 5*x^93 + x^94 + 2//7*x^95 + 5//2*x^97 + 8*x^99 + O(x^100))*y^2 + (7//10 + 3//7*x^53 + 3//4*x^54 + 1//2*x^56 + 3*x^57 + x^58 + x^59 + 9//8*x^60 + 5*x^61 + 3//4*x^62 + 4//3*x^63 + 2*x^64 + 1//2*x^67 + x^68 + 2//9*x^69 + x^70 + 5*x^71 + 9//10*x^72 + 10*x^73 + 9//8*x^74 + 8//7*x^75 + 7//8*x^76 + 2*x^77 + 1//2*x^78 + 2//5*x^79 + 10//7*x^80 + 5*x^81 + 2//7*x^82 + 3//5*x^83 + 3//2*x^85 + 1//3*x^87 + 2*x^88 + 9//5*x^89 + 3//2*x^90 + 9//5*x^91 + 1//2*x^92 + 2*x^93 + x^94 + 3//10*x^95 + 7*x^96 + 9//2*x^98 + 10//7*x^99 + O(x^100))*y + 1 + 7*x^70 + x^72 + 8*x^73 + 9//7*x^74 + 3//2*x^75 + 2//5*x^76 + 5//9*x^77 + 8//3*x^78 + 2//3*x^79 + 5//4*x^80 + 2*x^81 + 5//7*x^82 + 4//3*x^83 + 7//4*x^84 + 2//5*x^85 + 2//3*x^86 + 5//3*x^87 + 5//6*x^88 + x^89 + x^90 + 9//7*x^91 + 5//3*x^92 + 9//2*x^93 + 6//5*x^94 + 1//3*x^95 + 7//8*x^96 + 3//2*x^97 + 7//3*x^98 + 6//5*x^99 + O(x^100)
B = y^100 + (7 + 1//8*x^31 + 10//7*x^32 + 2//9*x^33 + 9//7*x^35 + 1//9*x^36 + 3//5*x^37 + x^38 + 7//10*x^39 + 8*x^40 + 3//5*x^41 + x^42 + 10*x^43 + 5//8*x^44 + 7//4*x^45 + 1//2*x^46 + 3//2*x^47 + 7//9*x^48 + 3//4*x^49 + 3//8*x^50 + 8*x^51 + 4//5*x^52 + 4//5*x^53 + 1//3*x^54 + 1//6*x^56 + 1//2*x^57 + x^58 + 3*x^59 + 3//5*x^60 + 2*x^61 + 1//2*x^62 + 8//9*x^63 + 7//4*x^64 + 1//2*x^65 + 9//5*x^66 + 5//2*x^67 + 5//6*x^68 + 1//10*x^69 + 8*x^70 + 8//7*x^72 + 5*x^74 + 7*x^75 + 4//5*x^76 + 5//8*x^77 + x^78 + 7//5*x^80 + 9*x^82 + 2*x^83 + 5//3*x^84 + 2*x^85 + 5//9*x^86 + 7*x^88 + 1//4*x^89 + 1//8*x^90 + x^91 + 2//9*x^92 + 2*x^93 + 9*x^94 + 3//2*x^95 + 4//3*x^96 + 3*x^97 + 2*x^98 + 1//7*x^99 + O(x^100))*y^12 + (1//8 + 4//3*x^7 + 5//2*x^8 + 4//9*x^9 + 5//2*x^10 + 2*x^11 + 5//6*x^13 + 4//9*x^14 + 1//2*x^15 + 1//6*x^16 + 7*x^17 + 7//4*x^18 + 8//3*x^19 + 6*x^20 + 2*x^22 + 10*x^23 + x^24 + 2*x^25 + 1//8*x^26 + 2//3*x^27 + 5//4*x^28 + 9//5*x^29 + 4//9*x^30 + 1//4*x^31 + x^32 + 7//4*x^33 + 10//3*x^34 + 7//10*x^35 + 4//5*x^36 + 3*x^37 + 1//8*x^38 + 4//3*x^39 + 3//2*x^40 + x^41 + 1//4*x^42 + x^43 + x^44 + 10//7*x^46 + 7//6*x^47 + x^48 + 4//9*x^50 + 4//5*x^51 + x^53 + 3//2*x^54 + 1//2*x^55 + x^57 + 4//9*x^58 + x^59 + 8//7*x^60 + x^62 + x^63 + 3//2*x^64 + x^66 + 5*x^67 + 1//6*x^69 + 9//7*x^70 + 2*x^71 + 3//2*x^72 + 7//10*x^73 + 10//3*x^74 + 7//5*x^75 + 4//3*x^76 + 1//9*x^77 + 5//9*x^78 + 5//6*x^79 + 3//8*x^80 + 7//9*x^81 + 2//3*x^82 + 8//3*x^83 + 3//4*x^84 + 1//8*x^85 + 7*x^87 + 8//7*x^88 + 5//7*x^89 + 7//10*x^90 + 1//3*x^91 + 7//9*x^92 + 7//9*x^93 + 6*x^94 + 2//3*x^95 + 1//3*x^96 + 3//5*x^97 + 7//2*x^98 + 1//2*x^99 + O(x^100))*y^11 + (1//2 + 4*x^51 + 1//3*x^53 + 7//8*x^54 + 7//3*x^55 + x^56 + 8//9*x^57 + 6//5*x^58 + 2*x^59 + x^60 + 6//5*x^61 + 4//5*x^62 + 5*x^63 + 8//7*x^65 + 1//2*x^66 + 3//7*x^67 + 9//4*x^68 + 5//4*x^69 + 1//9*x^70 + 1//8*x^71 + 5//2*x^72 + x^73 + 5//7*x^74 + 1//3*x^75 + 9//10*x^76 + 2//3*x^77 + 1//2*x^78 + 3*x^79 + 6//5*x^80 + x^81 + 5*x^82 + 6//7*x^84 + 3//2*x^85 + 5//4*x^86 + x^88 + 10//9*x^89 + 1//9*x^90 + 5//9*x^91 + 1//2*x^92 + x^93 + 3//5*x^94 + 3//10*x^95 + 9//8*x^96 + 1//3*x^99 + O(x^100))*y^10 + (7//10 + 2//5*x^75 + 1//2*x^76 + 1//3*x^77 + 2*x^78 + 1//5*x^79 + 6*x^80 + 7//2*x^81 + 5//7*x^82 + 7//5*x^83 + 4*x^84 + x^85 + 5//4*x^86 + x^87 + 4//3*x^88 + 5//4*x^89 + 5//9*x^90 + 9//7*x^91 + 3//2*x^92 + 9//10*x^93 + 5//4*x^94 + 10//7*x^95 + 5*x^96 + x^97 + 1//8*x^98 + 7//10*x^99 + O(x^100))*y^9 + (4//3 + 3//2*x^55 + 7//4*x^56 + 5//4*x^57 + x^58 + 3//7*x^59 + 1//5*x^60 + 4//5*x^61 + 8//5*x^62 + 8//7*x^64 + 4*x^65 + x^66 + 6*x^67 + x^68 + 7//5*x^69 + 9//4*x^70 + 6//7*x^71 + 9//8*x^72 + 10//9*x^73 + 9//10*x^74 + x^75 + 2//3*x^76 + 2//5*x^77 + 2//5*x^78 + 8//7*x^79 + 1//2*x^80 + 7//4*x^81 + 3//7*x^82 + x^83 + 8*x^84 + 2//5*x^85 + 6//5*x^86 + x^87 + 1//10*x^88 + 5*x^89 + 10//3*x^91 + 6//7*x^92 + 8//3*x^93 + 1//3*x^95 + 1//5*x^96 + 1//4*x^97 + 1//2*x^98 + 1//3*x^99 + O(x^100))*y^8 + (2 + 3//2*x^90 + 1//2*x^91 + 1//2*x^92 + 8//5*x^93 + 2//3*x^94 + 7//9*x^95 + 9//5*x^96 + 9//8*x^97 + 5//8*x^99 + O(x^100))*y^7 + (1 + 2*x^44 + 3//2*x^45 + 8//9*x^46 + 7//3*x^47 + 1//2*x^48 + 1//2*x^49 + 7*x^50 + 1//2*x^51 + 7//6*x^52 + 1//2*x^53 + 4//7*x^54 + 4*x^55 + 7//9*x^56 + 1//5*x^58 + 8*x^59 + 2//5*x^61 + 9//5*x^62 + 2*x^63 + 2//3*x^64 + 7//2*x^65 + 8//7*x^66 + 7//5*x^67 + 1//2*x^69 + 7//9*x^70 + 2*x^71 + 1//7*x^72 + 4//3*x^74 + 2//5*x^75 + 5//7*x^76 + 7//4*x^77 + 3//2*x^78 + 8//5*x^79 + 6//7*x^80 + 2//3*x^81 + 5//8*x^82 + 5//2*x^83 + 3//4*x^84 + 5//4*x^85 + 5//7*x^86 + 9//7*x^87 + 3//8*x^88 + 3//10*x^89 + 1//10*x^90 + x^91 + 7//4*x^92 + 2*x^93 + 9*x^94 + 2*x^95 + 10*x^96 + 1//2*x^97 + x^98 + 4//7*x^99 + O(x^100))*y^6 + (5//8 + x^84 + x^85 + 1//4*x^86 + 3*x^87 + 7//5*x^88 + 8//5*x^90 + 1//5*x^91 + 7//2*x^92 + 2*x^93 + 4//7*x^96 + 2//5*x^97 + 1//2*x^98 + 1//7*x^99 + O(x^100))*y^5 + (2//5 + 8*x^76 + 5//9*x^77 + x^78 + 5//2*x^79 + 5//2*x^80 + 1//8*x^81 + 3//2*x^82 + 4//9*x^84 + 1//5*x^85 + 2*x^86 + 3*x^87 + 3//5*x^88 + 5//7*x^89 + 8//7*x^90 + 5*x^91 + 7//4*x^92 + x^93 + 8*x^94 + 5*x^95 + 7//6*x^97 + 2//3*x^98 + x^99 + O(x^100))*y^4 + (5//4 + x^66 + 3//4*x^67 + x^70 + 5//4*x^71 + 5//7*x^72 + 6//5*x^73 + 2//3*x^75 + 9//7*x^76 + x^77 + 10//9*x^78 + 3//2*x^79 + x^80 + 1//3*x^82 + 8//7*x^83 + x^84 + 1//6*x^85 + 6//7*x^86 + 5*x^87 + 5//2*x^88 + 8//7*x^89 + 2//3*x^91 + 1//2*x^92 + 7//5*x^93 + 1//8*x^94 + 1//8*x^95 + 5*x^96 + 7//8*x^98 + 8//3*x^99 + O(x^100))*y^3 + (2//5 + 10*x^33 + 10//7*x^35 + 1//8*x^39 + 1//7*x^40 + 3//7*x^41 + x^42 + 5//3*x^43 + 3//2*x^44 + 1//4*x^45 + 5//8*x^46 + 5//4*x^47 + 2//3*x^48 + 9//10*x^49 + x^50 + 2//9*x^51 + 1//9*x^52 + 2//3*x^53 + 9//2*x^54 + 7//3*x^55 + 4//5*x^56 + 1//6*x^57 + 5//2*x^58 + 7//10*x^59 + 3//5*x^60 + 4//3*x^61 + x^62 + 10//9*x^63 + 1//10*x^64 + 4//5*x^66 + 4//3*x^67 + 9//5*x^68 + 5//6*x^70 + 1//3*x^71 + 5//2*x^72 + 5//4*x^73 + 1//5*x^74 + 3//4*x^75 + 5//9*x^76 + 5//2*x^77 + 10//9*x^78 + 3//10*x^80 + 9*x^81 + 1//9*x^82 + 1//3*x^83 + 10//9*x^84 + 9//2*x^85 + 6//7*x^86 + 4//3*x^87 + 7//2*x^88 + 1//3*x^89 + 2//3*x^90 + 1//9*x^91 + 5//2*x^92 + 3//7*x^93 + 1//4*x^94 + 1//5*x^95 + 2*x^98 + 8//5*x^99 + O(x^100))*y^2 + (5//3 + 9*x^23 + 8//3*x^24 + 5//9*x^25 + 3//10*x^26 + 10//7*x^27 + x^28 + 5//4*x^29 + 1//7*x^30 + 1//5*x^31 + 9*x^32 + 3//2*x^34 + x^35 + 10//7*x^36 + x^37 + 3//2*x^38 + 1//5*x^39 + 1//2*x^40 + 2//3*x^41 + 2*x^42 + x^43 + 3//4*x^44 + 5//3*x^45 + 10//9*x^46 + 9//4*x^47 + 1//3*x^48 + 2*x^49 + 1//4*x^50 + 1//2*x^51 + 6//5*x^52 + 9//4*x^53 + 5//3*x^54 + 7//3*x^56 + 7//6*x^58 + 1//2*x^59 + 5//6*x^60 + 3//5*x^61 + 5//4*x^62 + 3//5*x^63 + 1//2*x^64 + x^65 + 2//5*x^66 + 2*x^67 + x^68 + 5*x^69 + x^70 + 2*x^71 + 7//9*x^72 + 4//5*x^73 + x^75 + 4//3*x^76 + 2*x^77 + 5*x^78 + 4//3*x^79 + 3//2*x^81 + 7//2*x^82 + 8*x^83 + 2//3*x^84 + 5//9*x^85 + 9//5*x^86 + 10*x^87 + 1//2*x^88 + 5//3*x^89 + 1//3*x^90 + 4//3*x^91 + 2//9*x^92 + x^94 + 2*x^95 + 9//10*x^96 + 3//5*x^97 + x^98 + 5//3*x^99 + O(x^100))*y + 4//3 + 5//4*x^64 + 2*x^65 + 3//10*x^66 + 7//5*x^67 + 3//5*x^68 + 7//9*x^69 + 2//3*x^70 + 7//5*x^71 + x^72 + 4//9*x^73 + 3//2*x^74 + 8//7*x^76 + 1//2*x^77 + 4//5*x^78 + 2*x^79 + 3*x^80 + x^81 + 2*x^83 + 8//3*x^85 + 7//6*x^86 + 3//4*x^87 + 4//7*x^88 + 3*x^89 + 2//5*x^90 + 10//3*x^91 + 1//7*x^92 + 2*x^93 + 8//7*x^94 + 1//2*x^95 + 2//7*x^96 + 2//7*x^97 + 8//7*x^98 + 6//7*x^99 + O(x^100)

(Q,R) = fast_div_rem(A,B, l)

RX, _ = power_series_ring(RDF, l, "x"; model=:capped_absolute)
RXY, _ = polynomial_ring(RX,"y")

AF = to_other_poly(A, RDF, RX, RXY)
BF = to_other_poly(B, RDF, RX, RXY)

(QF, RF, r_Q, r_R) = val_fast_div_rem(AF, BF, l, rho, tau)

Q_prime = to_other_poly(Q, RDF, RX, RXY)
actual_Q_bound = biv_norm(Q_prime - QF, rho, tau)
given_Q_bound = mag(r_Q)

R_prime = to_other_poly(R, RDF, RX, RXY)
actual_R_bound = biv_norm(R_prime - RF, rho, tau)
given_R_bound = mag(r_R)

F = AF
G = BF

K = parent(G)
m = F.length - 1
n = G.length - 1

if m < n
    return (K(0), F, 0, 0)
end

FB = to_ball(F, l)
GB = to_ball(G, l)

(H_tilda, r_h) = val_fast_inv(rev(G, n), l, m-n+1, rho, tau)
H_tildaB = to_ball(H_tilda, l)

Q_tildaB = rev(biv_mullow(H_tildaB, rev(FB, m), l, m-n+1),m-n)
Q_tilda = from_ball(Q_tildaB, l)
error_Q = biv_norm(Q_tildaB - to_ball(Q_tilda, l), rho, tau)
r_Q = biv_norm(FB, rho, tau)*r_h + error_Q

FB
biv_truncate(FB, l, n)
R_tildaB = biv_truncate(FB, l, n) - biv_mullow(Q_tildaB, GB, l, n)
R_tilda = from_ball(biv_truncate(FB, l, n) - biv_mullow(Q_tildaB, GB, l, n), l)
error_R = biv_norm(R_tildaB - to_ball(R_tilda, l), rho, tau)
r_R = biv_norm(GB, rho, tau)*r_Q + error_R
RF


##################################################################################

QX, x = power_series_ring(AbstractAlgebra.QQ, 10, "x"; model=:capped_absolute)
QXY, y = polynomial_ring(QX,"y")
l = 10

g = y^15 + (4*x^9 + 5*x^3)*y^2 + x^6
h = y^15 + (324x^4 + 234*x+5) + 17* y^7
F = g*h


h = biv_truncate(h, 1)


K = parent(F)
G = biv_truncate(g, 1)
H = biv_truncate(h, 1)
S,T = cofactors(g, h, AbstractAlgebra.QQ, QX, QXY)
dG = G.length

S*G + T*H

# 
# Caste à cause de l'évaluation en 0
Hinv = K(rev(H, H.length-1)(0))
Q = K(0)

i = 1


i *= 2
    
if i > l
    i = l
end

# Update E
E = biv_truncate(F, i) - biv_mullow(G, H, i)

# Compute Q and R iteratively
SE = biv_mullow(S, E, i)
m = max(SE.length - 1, 0)
n = H.length - 1
m-n+1
# fast_div_rem(SE, H, i)
# SE = Q*H + R
if m >= n
    # Calcul de Q et R pour mettre à jour G et H
    Q = rev(biv_mullow(Hinv, rev(SE, m), i, m-n+1), m-n)
end
R = biv_truncate(SE, i, n) - biv_mullow(Q, H, i, n) 

# MAJ G et H
G = biv_truncate(G, i, dG) + biv_mullow(T, E, i, dG) + biv_mullow(Q, G, i, dG)
H = biv_truncate(H + R, i)

# fast_inv(rev(H,n), i, m-n+1)
#
Hinv = biv_mullow(Hinv, K(2) - rev(H,n)*Hinv, i, m-n+1)



(Qe,Re) = fast_div_rem(biv_mullow(S, E, i), H, i)


# Update S and T
E2 = biv_mullow(S, G, i) + biv_mullow(T ,H, i) - K(1)
(Q2,R2) = fast_div_rem(biv_mullow(S, Ehensel_version, ma2, i), H, i)
S = biv_truncate(S-R2, i)
T = biv_truncate(T,i, dG) - biv_mullow(T, E2, i, dG) - biv_mullow(Q2, G, i, dG)




###################################################################################
l = 10
QX, x = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
QXY, y = polynomial_ring(QX,"y")

A = (7//6 + 4//5*x^40 + 4//7*x^41 + 3*x^42 + 7//5*x^43 + 7//2*x^44 + 3//10*x^45 + 8*x^46 + 1//2*x^47 + 9//2*x^48 + 5*x^49 + O(x^50))*y^16 + (10//9 + 7//4*x^18 + x^19 + 9//4*x^20 + 8//9*x^21 + 7//3*x^22 + 1//3*x^23 + 5//2*x^24 + 1//3*x^25 + 2*x^26 + 5*x^28 + 1//5*x^29 + 5//8*x^30 + 3*x^31 + 3//2*x^33 + 1//3*x^34 + 5//4*x^35 + 1//2*x^36 + 9//5*x^37 + 3*x^38 + 3//8*x^39 + 1//7*x^40 + 1//5*x^42 + 8//7*x^43 + 5//3*x^44 + 9//2*x^45 + 3//2*x^46 + 7//2*x^47 + 7//8*x^48 + 2//7*x^49 + O(x^50))*y^14 + (5//3 + 2*x^32 + 9//7*x^33 + 2//3*x^34 + 3//8*x^36 + 4*x^37 + 9//8*x^38 + 3//5*x^39 + 7//2*x^41 + 3//2*x^42 + 9//8*x^43 + 5//8*x^44 + 1//2*x^45 + 9//8*x^46 + 3//4*x^47 + 5*x^48 + 1//2*x^49 + O(x^50))*y^13 + (1 + 2//3*x^41 + 2//5*x^42 + 3//4*x^43 + 4//3*x^44 + 8//9*x^45 + 6*x^46 + 6//7*x^47 + 5//8*x^48 + x^49 + O(x^50))*y^12 + (3//7 + 1//9*x^3 + 10//3*x^4 + 3//5*x^5 + 7//9*x^6 + 1//9*x^7 + x^8 + 5//4*x^9 + 5//6*x^10 + 2//3*x^11 + 3//4*x^12 + 1//2*x^13 + x^14 + 7//3*x^15 + 5//3*x^16 + 3*x^17 + 1//7*x^19 + 6//5*x^20 + 9//5*x^21 + 5//3*x^22 + 3//4*x^23 + 7//4*x^24 + 3//2*x^25 + 6*x^26 + 2//3*x^27 + 5//2*x^29 + 9//7*x^31 + 9//2*x^32 + x^33 + 3//2*x^34 + 1//3*x^35 + 9*x^36 + 9//4*x^37 + 5//6*x^38 + 4//3*x^39 + 7*x^40 + 5//2*x^41 + 1//4*x^42 + 6*x^43 + 3//4*x^44 + x^45 + 4//7*x^46 + 8//9*x^47 + 7//2*x^48 + 2//3*x^49 + O(x^50))*y^11 + (1//2 + 1//9*x^19 + 8//9*x^20 + 9//4*x^21 + 2*x^22 + 9//5*x^23 + x^25 + 5//7*x^26 + 6//7*x^27 + 4//5*x^28 + 6*x^30 + 1//7*x^31 + 3//7*x^33 + 3//7*x^35 + 1//3*x^36 + 5//2*x^37 + 3//2*x^38 + x^40 + 1//3*x^41 + 1//5*x^42 + 2//5*x^43 + 5//8*x^44 + 6*x^45 + 1//6*x^46 + 3//5*x^47 + 7//4*x^49 + O(x^50))*y^10 + (3//7 + 9*x^23 + 10//7*x^24 + 3//2*x^25 + 7//4*x^26 + 3//8*x^27 + 1//2*x^28 + 5//2*x^29 + 3//5*x^30 + 9//4*x^31 + 3//2*x^35 + 5//6*x^36 + 2//3*x^37 + x^38 + 3//4*x^39 + 7//4*x^40 + 3//2*x^41 + 1//5*x^42 + 3*x^43 + 1//3*x^44 + 8//5*x^45 + 4//5*x^46 + 7//3*x^47 + 1//5*x^48 + 9//4*x^49 + O(x^50))*y^9 + (4 + 4//7*x^41 + 10//9*x^43 + 10*x^44 + 4//5*x^45 + 2//3*x^46 + 7//5*x^47 + 2*x^48 + 5*x^49 + O(x^50))*y^8 + (2 + 5*x^31 + x^33 + 1//2*x^34 + 7*x^36 + 7//8*x^37 + 9*x^38 + 7//9*x^39 + 7//4*x^40 + 3//5*x^41 + 5*x^42 + 1//2*x^43 + 1//2*x^44 + x^46 + 3//2*x^47 + 1//2*x^48 + 3//4*x^49 + O(x^50))*y^7 + (3//7 + 5*x^24 + 6//5*x^26 + 5//8*x^27 + x^28 + 3//4*x^29 + 3//2*x^30 + 1//2*x^31 + x^32 + 3//5*x^33 + 3//7*x^34 + 1//2*x^35 + 3//10*x^36 + x^39 + 10//7*x^40 + 5//6*x^41 + 3//4*x^42 + 1//2*x^43 + 9//7*x^44 + 4//3*x^45 + 1//7*x^46 + 2//3*x^47 + 3*x^48 + 3//8*x^49 + O(x^50))*y^6 + (1//6 + O(x^50))*y^4 + (1//5 + 9//4*x^47 + 4//7*x^49 + O(x^50))*y^3 + (5//3 + 1//2*x^23 + x^24 + 5//4*x^25 + 4*x^26 + 5*x^27 + 3//4*x^28 + 4//3*x^30 + 2//3*x^31 + 5//8*x^32 + x^33 + 1//3*x^34 + 5//4*x^35 + 10//9*x^36 + 3*x^37 + x^39 + 1//2*x^40 + x^41 + 3//4*x^42 + 5//4*x^43 + x^44 + 5*x^45 + 2*x^46 + 8//5*x^47 + 8//9*x^48 + 4//7*x^49 + O(x^50))*y^2 + (3//2 + 4//5*x^7 + 4//9*x^8 + 3*x^9 + 7*x^11 + 7//8*x^13 + 7//4*x^14 + x^15 + 2*x^16 + 3*x^17 + 3//2*x^18 + x^19 + 5//2*x^20 + 9//10*x^21 + 1//6*x^22 + x^23 + 1//2*x^24 + 9//10*x^25 + 3//2*x^26 + 8//5*x^27 + 5//3*x^28 + 5//7*x^29 + 7//6*x^30 + 4//3*x^31 + 8//3*x^32 + x^33 + 2//3*x^34 + 4*x^35 + 3//5*x^36 + x^37 + 1//4*x^39 + 1//2*x^40 + 1//3*x^41 + 7*x^42 + 6//7*x^43 + 7//8*x^45 + 2*x^46 + 1//6*x^47 + 7//2*x^48 + 9//4*x^49 + O(x^50))*y + 1//5 + 9//2*x^34 + 5//8*x^35 + 3//4*x^36 + 5//3*x^37 + 1//2*x^38 + 7//6*x^39 + 4//5*x^41 + 8*x^42 + 9//8*x^43 + x^44 + 10//7*x^45 + 2//3*x^46 + 5//2*x^47 + 1//5*x^48 + 3//2*x^49 + O(x^50)
B = (1//4 + 1//7*x + 2//3*x^2 + 2//3*x^3 + 1//8*x^4 + 2//5*x^6 + 3//4*x^7 + 1//4*x^8 + 1//2*x^9 + 3//2*x^10 + 3//4*x^11 + 2//7*x^13 + 1//2*x^14 + 7//10*x^15 + 1//6*x^17 + x^18 + 9//7*x^19 + 9//5*x^20 + 5//9*x^21 + 5*x^22 + 1//4*x^24 + 5*x^25 + 10*x^26 + 9*x^27 + 5//9*x^28 + 2*x^29 + 2//3*x^30 + 1//10*x^31 + 6*x^32 + 1//2*x^34 + 3*x^35 + 1//9*x^36 + 3//2*x^37 + 4//7*x^38 + 3//7*x^39 + x^43 + 1//8*x^44 + 2*x^45 + 5//9*x^46 + 4//3*x^47 + 5*x^48 + 6//5*x^49 + O(x^50))*y^3 + (1 + 3*x^32 + 3*x^33 + 7//6*x^34 + 3//2*x^35 + 1//2*x^36 + 2*x^37 + 4*x^38 + 5*x^39 + 10//9*x^40 + 9//10*x^41 + 10*x^42 + 10//3*x^43 + 2//3*x^46 + 3//4*x^47 + 9//10*x^48 + x^49 + O(x^50))*y^2 + (1//8 + x + 7//3*x^2 + 1//2*x^3 + 4//3*x^6 + 8//5*x^7 + 9//4*x^8 + 5//4*x^9 + 7//2*x^10 + 4*x^12 + 8//5*x^13 + x^14 + 4*x^15 + 1//5*x^16 + x^17 + 3*x^18 + 8//3*x^19 + 4//7*x^20 + 1//5*x^22 + 5//6*x^23 + 2//3*x^24 + 8//9*x^26 + 4//3*x^27 + x^28 + 7//4*x^29 + x^30 + 7//3*x^31 + 2*x^34 + 5*x^35 + 1//2*x^36 + 4//5*x^37 + 5//2*x^38 + 7//10*x^39 + 1//2*x^40 + 9//4*x^41 + 3//8*x^42 + x^43 + 8//7*x^44 + 9//8*x^47 + 4//3*x^48 + 1//8*x^49 + O(x^50))*y + 5 + 9//8*x^47 + 2//5*x^48 + 10//9*x^49 + O(x^50)

t = @elapsed fast_div_rem(A,B, 10)

t

(Q,R) = fast_div_rem(A,B,50)

B = random_polynomial(AbstractAlgebra.QQ, QX, QXY, A.length - 1, l; monic=true)



A = random_polynomial(AbstractAlgebra.QQ, QX, QXY, 10, l)

gi = random_polynomial(AbstractAlgebra.QQ, QX, QXY, 10, 10)