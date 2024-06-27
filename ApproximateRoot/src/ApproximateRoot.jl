module ApproximateRoot

using AbstractAlgebra

export approximate_root
export nthroot

function approximate_root(F::Generic.Poly{T}, N::Int, Version::Int=1) where {T}
    if N == 1
        return P
    else
        G = reverse(F)
        if Version == 1
            p,_ = nthroot(G, N)
            return reverse(p)
        else
            p,_ = _nthroot_2(G, N)
            return reverse(q)
        end
    end
end


function nthroot(G::Generic.Poly{T}, N::Int) where {T}
    d = G.length - 1
    if d % N != 0
        error("N should divide the degree of the input polynomial")
    end

    m = div(d,N)

    # We make two Newton iterations at the same time to directly find the approximate root of P.
    # First iteration is y_{k+1} = y_k - (y_k - Q*z_k)/N 
    # Second one is z_{k+1} = 
    K = parent(G)
    y = K(G(0))
    z = K(1/G(0))
    k=1
    while k < m+1
        k *= 2
        y -= (y - mullow(G, z, k))/N
        z = truncate(2*z - y^(N-1)*z*z, k)
    end

    return truncate(y, m+1), truncate(z, m+1)
    
end

function _nthroot_2(G::Generic.Poly{T}, N::Int) where {T}
    d = G.length - 1
    if d % N != 0
        error("N should divide the degree of the input polynomial")
    end

    m = div(d,N)
    K = parent(G)

    # We first make a Newton iteration to find the inverse of the N-th approximate root
    # The recursive formula is y_{k+1} = y_k + (y_k - P*y_k^N+1)/N
    y = K(1 / Q(0))
    k=1
    while k < m+1
        k *= 2
        y += (y - mullow(P, y^(N+1), k))/N
    end

    # We then make a second Newton iteration to fin the inverse of y ie. the approximate root of Q
    # The recursive formula is z_{k+1} = 
    z = K(y(0)) # = Q(0)
    l=1
    while l < m + 1
        l *= 2
        z = truncate(2*z - y*z*z, l)
    end

    return truncate(z, m+1), truncate(y, m+1)

end


end # module ApproxRoot