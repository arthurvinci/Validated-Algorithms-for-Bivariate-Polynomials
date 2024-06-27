using AbstractAlgebra
using ValidatedHensel



g = y^15 + (4*x^9 + 5*x^3)*y^2 + x^6
h = y^15 + (324x^4 + 234*x+5) + 17* y^7
F = g*h
g = biv_truncate(g, 1)
h = biv_truncate(h, 1)
s,t = cofactors(g, h, AbstractAlgebra.QQ, QX, QXY)
s*g + t*h


G, H, S,T = fast_hensel_lifting(F,g,h,s,t,l)
K = parent(F)
G = K(g)
H = K(h)
S = K(s)
T = K(t)
m = G.length - 1
n = H.length - 1

# Raise Hinv to correct y precision
Hinv = fast_inv(rev(H, n), 1, m+n+1)
biv_truncate(Hinv*rev(H,n), 1, m+n+1)

i=1

    i *= 2

    if i > l
        i = l
    end
   

    # Factorization error
    E = biv_truncate(F, i) - biv_mullow(G, H, i)
    
    # Update Q and R
    SE = biv_mullow(S, E, i)
    d = degree(SE)

    Q = rev(biv_mullow(Hinv, rev(SE, d), i, d-n+1), d-n)
    R = biv_truncate(SE, i, n) - biv_mullow(Q, H, i, n)


    Q,R
    qt,rt = fast_div_rem(SE, H, i) 

    biv_truncate(SE - Q*H - R, i)
    biv_truncate(SE - qt*H -rt, i)
    # Update G, H, Hinv (because we use updated H and Hinv later)
    G += biv_mullow(T, E, i, m) + biv_mullow(Q, G, i, m)
    H += R
    Hinv = biv_mullow(Hinv, K(2) - biv_mullow(Hinv, rev(H,n), i, m+n+1), i, m+n+1)
    biv_truncate(Hinv*rev(H,n), 1, m)
    biv_truncate(F - G*H, i)

    # Cofactord error
    E2 = biv_mullow(S, G, i) + biv_mullow(T ,H, i) - K(1)

    # Update Q' and R'
    SE2 = biv_mullow(S, E2, i)
    d2 = degree(SE2)
    degree(H)
    Q2 = rev(biv_mullow(Hinv, rev(SE2, d2), i, d2-n+1), d2-n)
    R2 = biv_truncate(SE2, i, n) - biv_mullow(Q2, H, i, n)

    Q2,R2
    qt2,rt2 = fast_div_rem(SE2, H, i) 

    biv_truncate(SE2 - Q2*H - R2, i)
    biv_truncate(SE2 - qt2*H -rt2, i)

    # Update S and T
    S -= R2
    T -= biv_mullow(T, E2, i, m) + biv_mullow(Q2, G, i, m)
    biv_truncate(S*G + T*H - K(1), i)
    





biv_truncate(F -g*h, i, j) 
biv_truncate(s*g + t*h - 1, i, j)

gx = biv_truncate(g,1)
hx = biv_truncate(h,1)
sx = biv_truncate(s,1)
tx = biv_truncate(t,1)






#########################################################
l = 20
QX, x = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
QXY, y = polynomial_ring(QX,"y")

g = x^4*y^15 + (4*x^9 + 5*x^5)*y^2 + x^8 + 1
h = y^15
F = g*h
g = biv_truncate(g, 1)
h = biv_truncate(h, 1)
s,t = cofactors(g, h, AbstractAlgebra.QQ, QX, QXY)

s*g + t*h

@time fast_hensel_lifting(F,g,h,s,t,l)
@time hensel_lifting(F,g,h,s,t,l)


K = parent(F)
G = K(g)
H = K(h)
S = K(s)
T = K(t)
dF = F.length
n = h.length - 1
m = dF - n - 1

# Raise Hinv to correct y precision
Hinv = fast_inv(rev(H, n), 1, dF+m)
i=1

    old_i = i
    i *= 2

    if i > l
        i = l
    end


    # Factorization error
    E = shift_coeffs_right(biv_truncate(F, i) - biv_mullow(G, H, i), old_i)

    # Update Q and R
    SE = biv_mullow(S, E, i)
    d = degree(SE)
    Q = rev(biv_mullow(Hinv, rev(SE, d), i, d-n+1), d-n)
    R = biv_truncate(SE, i, n) - biv_mullow(Q, H, i, n)

    (QE, RE) = fast_div_rem(SE, H, i)

    SE - QE*H - RE

    SE - Q*H - R

    # Update G, H, Hinv (because we use updated H and Hinv later)
    G += shift_coeffs_left(biv_mullow(T, E, i, m+1) + biv_mullow(Q, G, i, m+1), old_i)
    H += shift_coeffs_left(R, old_i)
    Hinv = biv_mullow(Hinv, K(2) - biv_mullow(Hinv, rev(H,n), i, dF), i, dF)

    # Cofactord error
    E2 = shift_coeffs_right(biv_mullow(S, G, i) + biv_mullow(T ,H, i) - K(1), old_i)

    # Update Q' and R'
    SE2 = biv_mullow(S, E2, i)
    d2 = degree(SE2)
    Q2 = rev(biv_mullow(Hinv, rev(SE2, d2), i, d2-n+1), d2-n)
    R2 = biv_truncate(SE2, i, n+1) - biv_mullow(Q2, H, i, n+1)

    (Q2E, R2E) = fast_div_rem(SE2, H, i)

    # Update S and T
    S -= shift_coeffs_left(R2, old_i)
    T -= shift_coeffs_left(biv_mullow(T, E2, i, m+1) + biv_mullow(Q2, G, i, m+1), old_i)
    
    biv_truncate(F - G*H, i)
    biv_truncate(S*G + T*H - K(1), i)

