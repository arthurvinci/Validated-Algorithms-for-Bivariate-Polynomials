export hensel_lifting, intermediate_hensel_lifting, fast_hensel_lifting

"""
    hensel_lifting(F::Generic.Poly{TY}, g::Generic.Poly{TY}, h::Generic.Poly{TY}, s::Generic.Poly{TY}, t::Generic.Poly{TY}, l::Int) where {TY}

Computes the bivariate Hensel lifting of `g` and `h` such that `F = g*h` up to a given precision `l`.

# Arguments
- `F::Generic.Poly{TY}`: The polynomial to compute the factors of.
- `g::Generic.Poly{TY}`: The first factor of `F` to lift.
- `h::Generic.Poly{TY}`: The second factor of `F` to lift.
- `s::Generic.Poly{TY}`: The cofactor of `g`.
- `t::Generic.Poly{TY}`: The cofactor of `h`
- `l::Int`: The precision level in x.

# Returns
- `(G, H, S, T)::Tuple{Generic.Poly{TY}, Generic.Poly{TY}, Generic.Poly{TY}, Generic.Poly{TY}}`: The refined factorization of `F` and the cofactors.

"""
function hensel_lifting(F::Generic.Poly{TY}, g::Generic.Poly{TY}, h::Generic.Poly{TY}, s::Generic.Poly{TY}, t::Generic.Poly{TY}, l::Int) where {TY}
    K = parent(F)
    G = K(g)
    H = K(h)
    S = K(s)
    T = K(t)
    m = F.length - h.length

    
    i=1
    while i<l
        i *= 2
        
        if i > l
            i = l
        end

        # Update G and H
        E = biv_truncate(F, i) - biv_mullow(G, H, i)
        (Q,R) = fast_div_rem(biv_mullow(S, E, i), H, i)
        G += + biv_mullow(T, E, i, m+1) + biv_mullow(Q, G, i, m+1)
        H += biv_truncate(R, i)

        # Update S and T
        E2 = biv_mullow(S, G, i) + biv_mullow(T ,H, i) - K(1)
        (Q2,R2) = fast_div_rem(biv_mullow(S, E2, i), H, i)
        S -= biv_truncate(R2, i)
        T -=  biv_mullow(T, E2, i, m+1) + biv_mullow(Q2, G, i, m+1)
    end
    
    return (G, H, S, T)
end

function intermediate_hensel_lifting(F::Generic.Poly{TY}, g::Generic.Poly{TY}, h::Generic.Poly{TY}, s::Generic.Poly{TY}, t::Generic.Poly{TY}, l::Int) where {TY}
    K = parent(F)
    G = K(g)
    H = K(h)
    S = K(s)
    T = K(t)
    m = F.length - h.length

    
    i=1
    while i<l
        
        if i > l
            i = l
        end

        # Update G and H
        E =  shift_coeffs_right(biv_truncate(F, 2*i) - biv_mullow(G, H, 2*i), i)
        (Q,R) = fast_div_rem(biv_mullow(S, E, i), H, i)
        G += + shift_coeffs_left(biv_mullow(T, E, i, m+1) + biv_mullow(Q, G, i, m+1), i)
        H += shift_coeffs_left(R, i)

        # Update S and T
        E2 =  shift_coeffs_right(biv_mullow(S, G, 2*i) + biv_mullow(T ,H, 2*i) - K(1), i)
        (Q2,R2) = fast_div_rem(biv_mullow(S, E2, i), H, i)
        S -= shift_coeffs_left(R2, i)
        T -=  shift_coeffs_left(biv_mullow(T, E2, i, m+1) + biv_mullow(Q2, G, i, m+1), i)
        
        i *= 2
    end
    
    return (G, H, S, T)
end

function fast_hensel_lifting(F::Generic.Poly{TY}, g::Generic.Poly{TY}, h::Generic.Poly{TY}, s::Generic.Poly{TY}, t::Generic.Poly{TY}, l::Int) where {TY}
    K = parent(F)
    G = K(g)
    H = K(h)
    S = K(s)
    T = K(t)
    dF = F.length
    n = h.length - 1
    m = dF - n - 1

    # Raise Hinv to correct y precision
    Hinv = fast_inv(rev(H, n), 1, dF)
    i=1
    while i<l

        if i > l
            i = l
        end
    
        # Factorization error
        E = shift_coeffs_right(biv_truncate(F, 2*i) - biv_mullow(G, H, 2*i), i)
        
        # Update Q and R
        SE = biv_mullow(S, E, i)
        d = degree(SE)
        Q = rev(biv_mullow(Hinv, rev(SE, d), i, d-n+1), d-n)
        R = biv_truncate(SE, i, n) - biv_mullow(Q, H, i, n)

        # Update G, H, Hinv (because we use updated H and Hinv later)
        G += shift_coeffs_left(biv_mullow(T, E, i, m+1) + biv_mullow(Q, G, i, m+1), i)
        H += shift_coeffs_left(R, i)
        Hinv = biv_mullow(Hinv, K(2) - biv_mullow(Hinv, rev(H,n), 2*i, dF+m), 2*i, dF+m)

        # Cofactord error
        E2 = shift_coeffs_right(biv_mullow(S, G, 2*i) + biv_mullow(T ,H, 2*i) - K(1), i)

        # Update Q' and R'
        SE2 = biv_mullow(S, E2, i)
        d2 = degree(SE2)
        Q2 = rev(biv_mullow(Hinv, rev(SE2, d2), i, d2-n+1), d2-n)
        R2 = biv_truncate(SE2, i, n+1) - biv_mullow(Q2, H, i, n+1)

        # Update S and T
        S -= shift_coeffs_left(R2, i)
        T -= shift_coeffs_left(biv_mullow(T, E2, i, m+1) + biv_mullow(Q2, G, i, m+1), i)    

        i *= 2
    end
    
    return (G, H, S, T)
 end

