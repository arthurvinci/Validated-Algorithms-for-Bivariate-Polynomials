export fast_inv, fast_div_rem

"""
    fast_inv(G::Generic.Poly{T}, k::Int=-1) where {T}

Computes the inverse of an element `G` of `C[[X]][Y]_n` modulo `k`.

# Arguments
- `G::Generic.Poly{T}`: The polynomial to calculate the inverse of.
- l::Int
- k::Int`: The modulo value (default: -1, meaning the degree of G).

# Returns
- `Generic.Poly{T}`: The inverse of `G` modulo `x^l, y^k`.
"""
function fast_inv(G::Generic.Poly{T}, l::Int, k::Int=-1) where {T}
    K = parent(G)
    
    if k == -1
        k = G.length - 1
    end

    H = K(1)
    i = 1
    while i < k
        i *= 2
        H = biv_mullow(H, K(2) - G*H, l, i)
    end

    return truncate(H, k)
end

"""
    fast_div_rem(F::Generic.Poly{T}, G::Generic.Poly{T}) where {T}

Calculate the quotient and remainder of the polynomial division `F/G`, where F and G are in `C[[X]][Y]_n`.

# Arguments
- `F::Generic.Poly{T}`: The dividend polynomial.
- `G::Generic.Poly{T}`: The divisor polynomial.

# Returns
- `(Generic.Poly{T}, Generic.Poly{T})`: A tuple containing the quotient and remainder.
"""
function fast_div_rem(F::Generic.Poly{T}, G::Generic.Poly{T}, l::Int) where {T}
    K = parent(F)
    m = F.length - 1
    n = G.length - 1

    if m < n
        return (K(0), F)
    end
    
    H = fast_inv(rev(G, n), l, m-n+1)
    Q = rev(biv_mullow(H, rev(F, m), l, m-n+1), m-n)
    R = biv_truncate(F, l, n) - biv_mullow(Q, G, l, n)    

    return (Q,R)
end