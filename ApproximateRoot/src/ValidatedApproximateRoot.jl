module ApproximateRoot

using AbstractAlgebra
using ApproximateRoot

export val_approx_root

CC = ComplexF64
CX, x = power_series_ring(CC, 10, "x")
CXY, y = polynomial_ring(CX,"y")

function val_approx_root(F::CXY, N::Int) where {T}
    return reverse(_val_nthroot(reverse(G), N))
end

function _val_nthroot(G::Generic.Poly{T}, N::Int) where {T}
    p_tilda, q_tilda = nthroot(G, N)
    delta = mullow(mullow(p_tilda, q_tilda, N), p_tilda^N - G) / N 
    lambda(r)
end

function _inf_norm(P::CXY) 
    max_ret = _inf_norm(coeff(P, 0))
    for i in 1:(P.length - 1)
        max_ret = max(max_ret, _inf_norm(coeff(P, i)))
    end
    return max_ret
end

function _inf_norm(P::CX)
    max_ret  = coeff(P, 0)
    for i in 1:(P.length -1)
        max_ret = max(max_ret, coeff(P, i))
    end
    return max_ret
end

function _inf_norm(P::CC)
    return abs(P)
end

function _cone_norm(P::CXY)
end

function _cone_norm(P::CX)
end

function _cone_norm(P: CC)
end
    

end