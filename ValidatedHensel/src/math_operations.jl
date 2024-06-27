export biv_mullow, biv_truncate, rev, to_ball, biv_norm, mag, from_ball, random_polynomial, to_other_poly, max_diff, cofactors, min_interval, max_interval,shift_coeffs_right, shift_coeffs_left, to_polynomial

"""
    biv_mullow(P::Generic.Poly{T}, Q::Generic.Poly{T}, l::Int, m::Int=-1) where {T}

Performs bivariate polynomial multiplication truncated to l inner terms and m outer terms.

# Arguments
- `P`: First polynomial.
- `Q`: Second polynomial.
- `l`: Truncation length for the inner coefficients.
- `m`: Degree to truncate the result. If not specified, it doesn't truncate.

# Returns
The truncated product of `P` and `Q`.
"""
function biv_mullow(P::Generic.Poly{T}, Q::Generic.Poly{T}, l::Int, m::Int=-1) where {T}

    if m==-1
        m = P.length + Q.length - 1
    end

    result = mullow(P,Q, m)

    inner_truncate!(result,l)
    return result
end


"""
    biv_truncate(P::Generic.Poly{T}, l::Int, m::Int=-1) where {T}

Truncates a bivariate polynomial.

# Arguments
- `P`: The polynomial to truncate.
- `l`: The inner truncation length.
- `m`: The degree to truncate `P`. If `-1`, the polynomial is not truncated by degree.

# Returns
The doubly truncated polynomial.
"""
function biv_truncate(P::Generic.Poly{T}, l::Int, m::Int=-1) where {T}
    if(m != -1)
        result = truncate(P, m)
    else
        result = deepcopy(P)
    end

    inner_truncate!(result,l)
    return result
end

function inner_truncate!(P::Generic.Poly{T}, l::Int) where {T}
    for i in 0:(P.length - 1)
        p_i = coeff(P, i)
        prec = p_i.prec
        new_coeff = truncate(p_i, l)
        new_coeff.prec = prec
        set_coefficient!(P, i, new_coeff)
    end
end

"""
    rev(P::Generic.Poly{T}, degree::Int) where T

Reverses the coefficients of the polynomial `P` up to a specified degree ie. computes rev_d(P)

# Arguments
- `P::Generic.Poly{T}`: The polynomial whose coefficients will be reversed.
- `degree::Int`: The degree up to which the coefficients will be reversed.

# Returns
The polynomial with reversed coefficients up to the specified degree.
"""
function rev(P::Generic.Poly{T}, degree::Int) where {T}
    new_P = parent(P)(0)
    degree+=1
    fit!(new_P, degree)
    for i = 1:(degree)
        setcoeff!(new_P, i-1, coeff(P, degree -i))
    end
    return new_P 
end


"""
    to_ball(F)

Converts a bivariate polynomial to an interval polynomial.

# Arguments
- `F`: A bivariate polynomial.

# Returns
The ball arithmetic equivalent.
"""
function to_ball(F::Generic.Poly{Generic.AbsSeries{T}}, l::Int) where {T}
    CBX, _ = power_series_ring(CBall, l, "x"; model=:capped_absolute)
    CBXY, _ = polynomial_ring(CBX, "y")

    ball_poly = CBXY(0)
    f = x -> CBall(string(x))
    for i in 0:(F.length -1)
        new_coeffs = map_coefficients(f, coeff(F, i); parent=CBX)
        set_coefficient!(ball_poly, i, new_coeffs)
    end
    return ball_poly
end

"""
    to_ball(F::Generic.AbsSeries{T}, l::Int) where T

Converts a series `F` to a ball series up to a specified precision `l`.

# Arguments
- `F::Generic.AbsSeries{T}`: The series to be converted.
- `l::Int`: The precision up to which the coefficients will be truncated.

# Returns
The ball series obtained by converting the coefficients of `F` to `CBall` up to precision `l`.
"""
function to_ball(F::Generic.AbsSeries{T}, l::Int) where {T}
    CBX, _ = power_series_ring(CBall, l, "x"; model=:capped_absolute)

    f = x -> CBall(string(x))
    return map_coefficients(f, F; parent=CBX)
end

"""
    from_ball(F::Generic.Poly{Generic.AbsSeries{T}}, l::Int) where T

Converts a bivariate ball polynomial to a floating point bivariate polynomial up to a specified precision `l`.

Problème: On prend le centre de la boule 

# Arguments
- `F::Generic.Poly{Generic.AbsSeries{T}}`: The bivariate polynomial to be converted.
- `l::Int`: The precision up to which the coefficients will be truncated.

# Returns
The bivariate polynomial obtained by the converison up to precision `l`.
"""
function from_ball(F::Generic.Poly{Generic.AbsSeries{T}}, l::Int) where {T}
    CFX, _ = power_series_ring(CFloat, l, "x"; model=:capped_absolute)
    CFXY, _ = polynomial_ring(CFX, "y")

    float_poly = CFXY(0)
    f = x -> CFloatElem(midpoint(x))
    for i in 0:(F.length - 1)
        new_coeffs = map_coefficients(f, coeff(F,i); parent=CFX)
        set_coefficient!(float_poly, i, new_coeffs)
    end
    return float_poly
end

"""
    from_ball(F::Generic.AbsSeries{T}, l::Int) where T

Converts a ball series to a floating point series up to a specified precision `l`.

# Arguments
- `F::Generic.AbsSeries{T}`: The series to be converted.
- `l::Int`: The precision up to which the coefficients will be truncated.

# Returns
The series obtained by the converison up to precision `l`.
"""
function from_ball(F::Generic.AbsSeries{T}, l::Int) where {T}
    CFX, _ = power_series_ring(CFloat, l, "x"; model=:capped_absolute)
    f = x -> CFloatElem(midpoint(x))
    return map_coefficients(f, F; parent=CFX)
end

"""
    biv_norm(F::Generic.Poly{Generic.AbsSeries{T}}, rho::Float64, tau::Float64) where T

Computes the bivariate norm of the given bivariate polynomial for given convergence radiuses.

# Arguments
- `F::Generic.Poly{Generic.AbsSeries{T}}`: The bivariate polynomial.
- `rho::Float64`: convergence radius in `x`
- `tau::Float64`: optional convergence radius in `y` (1 by default).

# Returns
An interval with the computed norm.

"""
function biv_norm(F::Generic.Poly{Generic.AbsSeries{T}}, rho::Float64, tau::Float64=1.0) where {T}
    f = x -> biv_norm(x, rho)
    F_abs = map_coefficients(f, F)
    eval = evaluate(F_abs, tau)
    return CBall(Float64(abs(eval), RoundUp))
end


"""
    biv_norm(F::Generic.AbsSeries{T}, rho::Float64) where T

Computes the norm of the given series for the given convergence radius.

# Arguments
- `F::Generic.AbsSeries{T}`: The series.
- `rho::Float64`: convergence radius of the series.

# Returns
An interval with the computed norm.
"""
function biv_norm(F::Generic.AbsSeries{T}, rho::Float64) where {T}
    f = x -> abs(x)
    F_abs = map_coefficients(f, to_polynomial(F))
    return evaluate(F_abs, rho)
end

"""
    to_polynomial(F::Generic.AbsSeries{CBallElem})

Converts a series with coefficients in `CBallElem` to a polynomial with coefficients in `CBall`.

# Arguments
- `F::Generic.AbsSeries{CBallElem}`: The series with coefficients in `CBallElem`.

# Returns
The polynomial obtained by converting the coefficients of `F` to `CBall`.
"""
function to_polynomial(F::Generic.AbsSeries{CBallElem})
    F_poly = polynomial(CBall, [0])
    for i in 0:(F.length - 1)
        set_coefficient!(F_poly, i, coeff(F, i))
    end
    return F_poly
end

"""
    to_polynomial(F::Generic.AbsSeries{CFloatElem})

Converts a series with coefficients in `CFloatElem` to a polynomial with coefficients in `CFloat`.

# Arguments
- `F::Generic.AbsSeries{CFloatElem}`: The series with coefficients in `CFloatElem`.

# Returns
The polynomial obtained by converting the coefficients of `F` to `CFloat`.
"""
function to_polynomial(F::Generic.AbsSeries{CFloatElem})
    F_poly = polynomial(CFloat, [0])
    for i in 0:(F.length - 1)
        set_coefficient!(F_poly, i, coeff(F, i))
    end
    return F_poly
end

"""
    mag(x::T) where T <: Union{arb, acb}

Computes the magnitude of a given interval.

# Arguments
- `x::T`: The interval whose magnitude is to be computed.

# Returns
The magnitude of the interval.
"""
function mag(x::T) where T <: Union{arb,acb}
    return Float64(abs(x), RoundUp)
end

"""
    random_polynomial(BaseType::U, SeriesType::AbstractAlgebra.Generic.AbsPowerSeriesRing{T}, PolyType::AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.AbsSeries{T}}, y_deg::Int, x_deg::Int) where {T, U}

Generates a random polynomial with coefficients from specified types and degrees. 
Number coefficients are randomly generated between 0 and 10.

# Arguments
- `BaseType::U`: The base type of the polynomial coefficients.
- `SeriesType::AbstractAlgebra.Generic.AbsPowerSeriesRing{T}`: The type of power series for the coefficients.
- `PolyType::AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.AbsSeries{T}}`: The type of polynomial ring for the polynomial.
- `y_deg::Int`: The degree of the polynomial in the y variable.
- `x_deg::Int`: The degree of the series in the x variable.


# Returns
The randomly generated polynomial.

"""
function random_polynomial(BaseType::U, SeriesType::AbstractAlgebra.Generic.AbsPowerSeriesRing{T}, PolyType::AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.AbsSeries{T}}, y_deg::Int, x_deg::Int; monic::Bool=false) where {T, U}
    
    Kint, _ = polynomial_ring(BaseType, "y")

    poly = rand(Kint, 0:y_deg, 0:10)
    f = x -> SeriesType(x) + rand(SeriesType, 1:x_deg, 0:10)
    ret = map_coefficients(f, poly; parent=PolyType)

    if monic
        set_coefficient!(ret, y_deg, SeriesType(1))
    else
        leading_coeff = SeriesType(0)
            while leading_coeff == SeriesType(0)
                leading_coeff = rand(SeriesType, 1:x_deg, 0:10)
            end
            set_coefficient!(ret, y_deg, leading_coeff)
    end
end

"""
    to_other_poly(F::Generic.Poly{Generic.AbsSeries{U}}, BaseType::V, SeriesType::AbstractAlgebra.Generic.AbsPowerSeriesRing{T}, PolyType::AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.AbsSeries{T}}) where {T, U, V}

Converts a bivariate polynomial to another bivariate polynomial with coefficients of a different base type.

# Arguments
- `F::Generic.Poly{Generic.AbsSeries{U}}`: The polynomial to convert.
- `BaseType::V`: The base type of the coefficients to convert to.
- `SeriesType::AbstractAlgebra.Generic.AbsPowerSeriesRing{T}`: The type of power series for the coefficients.
- `PolyType::AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.AbsSeries{T}}`: The type of polynomial ring for the result.

# Returns
The polynomial obtained by converting the coefficients of `F` to another base type.

"""
function to_other_poly(F::Generic.Poly{Generic.AbsSeries{U}}, BaseType::V, SeriesType::AbstractAlgebra.Generic.AbsPowerSeriesRing{T}, PolyType::AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.AbsSeries{T}}) where {T, U, V}
    f = x -> to_other_poly(x, BaseType, SeriesType)
    return map_coefficients(f, F; parent=PolyType)
end

"""
    to_other_poly(F::Generic.AbsSeries{U}, BaseType::V, SeriesType::AbstractAlgebra.Generic.AbsPowerSeriesRing{T}) where {T, U, V}

Converts a serie to another series with coefficients of a different base type.

# Arguments
- `F::Generic.AbsSeries{U}`: The seriesto convert.
- `BaseType::V`: The base type of the coefficients to convert to.
- `SeriesType::AbstractAlgebra.Generic.AbsPowerSeriesRing{T}`: The type of power series for the coefficients.

# Returns
The series obtained by converting the coefficients of `F` to another base type.

"""
function to_other_poly(F::Generic.AbsSeries{U}, BaseType::V, SeriesType::AbstractAlgebra.Generic.AbsPowerSeriesRing{T}) where {T, U, V}
    f = x -> BaseType(x)
    return map_coefficients(f, F; parent=SeriesType)
end

"""
    cofactors(F::Generic.Poly{Generic.AbsSeries{U}}, G::Generic.Poly{Generic.AbsSeries{U}}, BaseType::V, SeriesType::AbstractAlgebra.Generic.AbsPowerSeriesRing{T}, PolyType::AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.AbsSeries{T}}) where {T, U, V}

Computes the cofactors of two bivariate polynomials by taking coefifcients modulo O(x^1).

# Arguments
- `F::Generic.Poly{Generic.AbsSeries{U}}`: The first bivariate polynomial.
- `G::Generic.Poly{Generic.AbsSeries{U}}`: The second bivariate polynomial.
- `BaseType::V`: The base type of the coefficients for the resulting polynomials.
- `SeriesType::AbstractAlgebra.Generic.AbsPowerSeriesRing{T}`: The type of power series for the coefficients.
- `PolyType::AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.AbsSeries{T}}`: The type of polynomial ring for the result.
# Returns
The cofactors of `F` and `G`.

"""
function cofactors(F::Generic.Poly{Generic.AbsSeries{U}}, G::Generic.Poly{Generic.AbsSeries{U}}, BaseType::V, SeriesType::AbstractAlgebra.Generic.AbsPowerSeriesRing{T}, PolyType::AbstractAlgebra.Generic.PolyRing{AbstractAlgebra.Generic.AbsSeries{T}}) where {T, U, V}
    new_K, _ = polynomial_ring(BaseType, "y")

    new_F = map_coefficients( x -> BaseType(coeff(x, 0)), F; parent=new_K)
    new_G = map_coefficients( x -> BaseType(coeff(x, 0)), G; parent=new_K)

    _, s, t = gcdx(new_F, new_G)

    new_s = map_coefficients(x -> SeriesType(x), s; parent=PolyType)
    new_t = map_coefficients(x -> SeriesType(x), t; parent=PolyType)
    
    return (new_s, new_t)
end

"""
    min_interval(x::CBallElem)

Computes the minimum of an interval .

# Arguments
- `x::CBallElem`: The interval.

# Returns
The minimum of the interval `x`.

"""
function min_interval(x::CBallElem)
    center = Float64(midpoint(x))
    rad = Float64(radius(x))
    return center - rad
end

function max_interval(x::CBallElem)
    center = Float64(midpoint(x))
    rad = Float64(radius(x))
    return center + rad
end

function shift_coeffs_right(F::Generic.Poly{Generic.AbsSeries{T}}, n::Int) where {T}
    return map_coefficients( x -> shift_right(x, n), F)
end

function shift_coeffs_left(F::Generic.Poly{Generic.AbsSeries{T}}, n::Int) where {T}
    return map_coefficients(x -> shift_left(x,n), F)
end