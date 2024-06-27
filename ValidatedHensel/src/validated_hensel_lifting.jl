export val_hensel_lifting

"""
    val_hensel_lifting(F::Generic.Poly{TY}, g::Generic.Poly{TY}, h::Generic.Poly{TY}, s::Generic.Poly{TY}, t::Generic.Poly{TY}, l::Int, rho::Float64, tau::Float64) where{TY}

Computes the validated bivariate Hensel lifting of `g` and `h` such that `F = g*h` up to a given precision `l`.

# Arguments
- `F::Generic.Poly{TY}`: The polynomial to compute the factors of.
- `g::Generic.Poly{TY}`: The first factor of `F` to lift.
- `h::Generic.Poly{TY}`: The second factor of `F` to lift.
- `s::Generic.Poly{TY}`: The cofactor of `g`.
- `t::Generic.Poly{TY}`: The cofactor of `h`
- `l::Int`: The precision level in x.
- `rho::Float64`: The convergence radius in `x`.
- `tau::Float64`: The convergence raidus in `y` (defaults to 1.0).

# Returns
- `(G, H, r)::Tuple{Generic.Poly{TY}, Generic.Poly{TY}, Float64}`: The refined factorization of `F` and the validated bounds.

"""
function val_hensel_lifting(F::Generic.Poly{TY}, g::Generic.Poly{TY}, h::Generic.Poly{TY}, s::Generic.Poly{TY}, t::Generic.Poly{TY}, l::Int, rho::Float64, tau::Float64) where{TY}
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
    norm_G = mag(biv_norm(GB, rho, tau))
    norm_H = mag(biv_norm(HB, rho, tau))
    norm_S = mag(biv_norm(SB, rho, tau))
    norm_T = mag(biv_norm(TB, rho, tau))
    

    EB = biv_truncate(FB, l) - biv_mullow(GB, HB, l)

    (QB, RB, r_Q, r_R) = val_fast_div_rem_interval(biv_mullow(SB, EB, l, m+n), HB, l, rho, tau)
    delta_G = mag(biv_norm(biv_truncate(TB*EB + QB*GB, l, m + 1), rho, tau)) + norm_G*r_Q
    delta_H = mag(biv_norm(RB, rho, tau)) + r_R
    
    delta = max(delta_G, delta_H)


    # Lipschitz ratio bounding
    Y_poly = K(0)
    set_coefficient!(Y_poly, m + 2*n - 2, k(1))

    (PHIB, _, r_PHI, _) = val_fast_div_rem_interval(to_ball(Y_poly, l), HB, l, rho, tau)
    
    PSIB = SB*GB + TB*HB - to_ball(K(1), l)
    
    mu = max(norm_G, norm_H)
    nu = max(norm_T, norm_S)
    phi = mag(biv_norm(PHIB, rho, tau)) + r_PHI
    psi = mag(biv_norm(PSIB, rho, tau))

    # Compute the polynomials' roots
    a = 2*nu*(1+mu*phi)
    b = psi*(1 + mu*phi) - 1
    c = delta
    disc = b*b - 4*a*c


    if is_negative(disc) # Check that the left side of the interval is positive
        error("Could not validate bounds (polynomial has no real roots)"* string(b)*string(a)*string(c))
    end

    r_min = (-b - sqrt(disc))/(2*a)
    r_max = (-b + sqrt(disc))/(2*a)

    r = r_min
    if !is_positive(r_min)
        if !is_positive(r_max)
            error("Could not validate bounds (polynomial has no positive roots)")
        end
        r = r_max
    end

    return (G, H, r)
end
