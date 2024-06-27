export val_fast_inv, val_fast_div_rem, val_fast_inv_interval, val_fast_div_rem_interval

"""
    val_fast_inv(G::Generic.Poly{Generic.AbsSeries{T}}, l::Int, k::Int, rho::Float64, tau::Float64=1.0) where {T}

Computes the validated inverse of a bivariate polynomial.

# Arguments
- `G::Generic.Poly{Generic.AbsSeries{T}}`: The bivariate polynomial to compute the inverse of.
- `l::Int`: The precision level in `x`.
- `k::Int`: Another precision level in `y`.
- `rho::Float64`: The convergence radius in `x`.
- `tau::Float64`: The convergence raidus in `y` (defaults to 1.0).

# Returns
- `(H_tilda, r)::Tuple{Generic.Poly{Generic.AbsSeries{T}}, Float64}`: The inverse of `G` and the validated bound `r`.

"""
function val_fast_inv(G::Generic.Poly{Generic.AbsSeries{T}}, l::Int, k::Int, rho::Float64, tau::Float64=1.0) where {T}
    K = parent(G)

    H_tilda = fast_inv(G, l, k)
    H_tildaB = to_ball(H_tilda, l)
    GB = to_ball(G, l)

    delta = biv_norm(H_tildaB - biv_mullow(GB, H_tildaB^2, l, k), rho, tau)
    lambda = biv_norm(to_ball(K(1), l) - biv_mullow(GB, H_tildaB, l, k), rho, tau)

    if max_interval(lambda) > 1
        error("lambda is greater than 1: "*string(lambda))
    end

    if !is_nonnegative(delta)
        error("delta is not positive: "*string(midpoint(delta))*" "*string(radius(delta)))
    end
    
    r = delta/(1-lambda)

    if isnan(mag(r))
        error("the fast inversion bound is a NaN")
    end

    if !is_nonnegative(r)
        error("the fast inversion bound is not positive ", r)
    end

    return (H_tilda, r)
end

function val_fast_inv_interval(GB::Generic.Poly{Generic.AbsSeries{T}}, l::Int, k::Int, rho::Float64, tau::Float64) where {T}
    K = parent(GB)

    G = from_ball(GB, l)

    H_tilda = fast_inv(G, l, k)
    H_tildaB = to_ball(H_tilda, l)

    delta = biv_norm(H_tildaB - biv_mullow(GB, H_tildaB^2, l, k), rho, tau)
    lambda = biv_norm(to_ball(K(1), l) - biv_mullow(GB, H_tildaB, l, k), rho, tau)

    if max_interval(lambda) > 1
        error("lambda is greater than 1: "*string(lambda))
    end

    if !is_nonnegative(delta)
        error("delta is not positive: "*string(delta))
    end

    r = delta/(1-lambda)

    if isnan(mag(r))
        error("the fast inversion bound is a NaN")
    end

    if !is_nonnegative(r)
        error("the fast inversion bound is not positive ", r)
    end

    return (H_tildaB, r)
end


"""
    val_fast_div_rem(F::Generic.Poly{Generic.AbsSeries{T}}, G::Generic.Poly{Generic.AbsSeries{T}}, l::Int, rho::Float64, tau::Float64) where {T}

Computes the validated euclidian division of two bivariate polynomials.

# Arguments
- `F::Generic.Poly{Generic.AbsSeries{T}}`: The dividend bivariate polynomial.
- `G::Generic.Poly{Generic.AbsSeries{T}}`: The divisor bivariate polynomial.
- `l::Int`: The precision level in `x`.
- `rho::Float64`: The convergence radius in `x`.
- `tau::Float64`: The convergence raidus in `y` (defaults to 1.0).

# Returns
The validated quotient and remainder of the euclidian division of `F` by `G` along with their respective validated bounds.

"""
function val_fast_div_rem(F::Generic.Poly{Generic.AbsSeries{T}}, G::Generic.Poly{Generic.AbsSeries{T}}, l::Int, rho::Float64, tau::Float64=1.0) where {T}
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
    
    R_tildaB = biv_truncate(FB, l, n) - biv_mullow(Q_tildaB, GB, l, n)
    R_tilda = from_ball(biv_truncate(FB, l, n) - biv_mullow(Q_tildaB, GB, l, n), l)
    error_R = biv_norm(R_tildaB - to_ball(R_tilda, l), rho, tau)
    r_R = biv_norm(GB, rho, tau)*r_Q + error_R

    if midpoint(r_Q) == 0
        r_R = radius(r_R)
    end

    if isnan(mag(r_R))
        error("the fast div rem R bound is a NaN")
    end

    if !is_nonnegative(r_R)
        error("the fast div rem R bound is not positive ", r_R)
    end


    if isnan(mag(r_Q))
        error("the fast div rem Q bound is a NaN")
    end

    if !is_nonnegative(r_Q)
        error("the fast div rem Q bound is not positive ", r_Q)
    end

    return(Q_tilda, R_tilda, r_Q, r_R)
end


function val_fast_div_rem_interval(FB::Generic.Poly{Generic.AbsSeries{T}}, GB::Generic.Poly{Generic.AbsSeries{T}}, l::Int, rho::Float64, tau::Float64) where {T}
    K = parent(GB)
    m = FB.length - 1
    n = GB.length - 1

    if m < n
        return (K(0), F, 0, 0)
    end

    (H_tildaB, r_h) = val_fast_inv_interval(rev(GB, n), l, m-n+1, rho, tau)
    
    Q_tildaB = rev(biv_mullow(H_tildaB, rev(FB, m), l, m-n+1),m-n)
    r_Q = biv_norm(FB, rho, tau)*r_h
    
    # Center can  equal 0 in certain cases and create issues. In this case, the error is the radius
    if midpoint(r_Q) == 0
        r_Q = radius(r_Q)
    end



    R_tildaB = biv_truncate(FB, l, n) - biv_mullow(Q_tildaB, GB, l, n)
    r_R = biv_norm(GB, rho, tau)*r_Q

    if midpoint(r_Q) == 0
        r_R = radius(r_R)
    end

    if isnan(mag(r_R))
        error("the fast div rem R bound is a NaN")
    end

    if !is_nonnegative(r_R)
        error("the fast div rem R bound is not positive ", r_R)
    end


    if isnan(mag(r_Q))
        error("the fast div rem Q bound is a NaN")
    end

    if !is_nonnegative(r_Q)
        error("the fast div rem Q bound is not positive ", r_Q)
    end

    
    return (Q_tildaB, R_tildaB, r_Q, r_R)
end