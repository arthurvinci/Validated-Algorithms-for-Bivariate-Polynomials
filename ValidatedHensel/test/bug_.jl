using AbstractAlgebra
using Nemo
using ValidatedHensel


l = 10

QX, x = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
QXY, y = polynomial_ring(QX,"y")
RX, _ = power_series_ring(RDF, l, "x"; model=:capped_absolute)
RXY, _ = polynomial_ring(RX,"y")

rho = 0.1
tau = 1.0
(Fe,ge,he,se,te) = random_hensel_setup(10, 10)
F = to_other_poly(Fe, RDF, RX, RXY)
g = to_other_poly(ge, RDF, RX, RXY)
h = to_other_poly(he, RDF, RX, RXY)
s = to_other_poly(se, RDF, RX, RXY)
t = to_other_poly(te, RDF, RX, RXY)


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
    norm_G = biv_norm(GB, rho, tau)
    norm_H = biv_norm(HB, rho, tau)
    norm_S = biv_norm(SB, rho, tau)
    norm_T = biv_norm(TB, rho, tau)
    

    EB = biv_truncate(FB, l) - biv_mullow(GB, HB, l)
    
    # (QB, RB, r_Q, r_R) = val_fast_div_rem_interval(biv_mullow(SB, EB, l, m+n), HB, l, rho, tau)
    
    FB = biv_mullow(SB, EB, l, m+n)
    GB = HB
    m = degree(FB)
    n = degree(GB)


    (H_tildaB, r_h) = val_fast_inv_interval(rev(GB, n), l, m-n+1, rho, tau)
    
    Q_tildaB = rev(biv_mullow(H_tildaB, rev(FB, m), l, m-n+1),m-n)
    r_Q = biv_norm(FB, rho, tau)
    r_h
    r_Q * r_h

    
    R_tildaB = biv_truncate(FB, l, n) - biv_mullow(Q_tildaB, GB, l, n)
    r_R = biv_norm(GB, rho, tau)*r_Q



    Ok = H_tildaB - biv_mullow(rev(GB,n), H_tildaB^2, l, m-n+1)
    Ok
    test = to_polynomial(coeff(Ok, 6))
    t2 = map_coefficients(abs, test)