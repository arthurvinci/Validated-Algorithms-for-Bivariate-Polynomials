using AbstractAlgebra
using Nemo
using ValidatedHensel
using GLM
using DataFrames
using DoubleFloats

DFloat = RDF

# Computes numeric solution found with fast hensel and rhos
function get_numeric_value_and_rhos(F::Generic.Poly{T}, g::Generic.Poly{T}, h::Generic.Poly{T}, s::Generic.Poly{T}, t::Generic.Poly{T}, l::Int) where {T}
    
    # Compute log search rho 
    (gv, hv, r, log_rho) = val_find_rho_log(F, g, h, s, t, l)

    # Compute the approximate rhos for gn and hn
    g_incr = get_exp_incr(gv)
    h_incr = get_exp_incr(hv)

    return (gv, hv, r, g_incr, h_incr, log_rho)
end


# Validates hensel lifting by logarithmic search on rho with 4 decimals of precision
function val_find_rho_log(F::Generic.Poly{T}, g::Generic.Poly{T}, h::Generic.Poly{T}, s::Generic.Poly{T}, t::Generic.Poly{T}, l::Int) where {T}
    max_rho = 1.0
    min_rho = 0
    gn = 0
    hn = 0
    r = 0
    iterations = 0
    while (max_rho-min_rho)/max_rho > 1.e-3 && iterations < 30
        new_rho = (max_rho+min_rho)/2
        try
            (gv, hv , rv) = val_hensel_lifting(F, g, h, s, t, l, new_rho, 1.0)
            gn = gv
            hn = hv
            r = rv
            min_rho = new_rho
        catch
            max_rho = new_rho
        end
        iterations += 1
    end

    if r == 0
        error("Could not find a working rho")
    end

    return (gn, hn, r, min_rho)
end

# Returns the approximate convergence radius by linear regression for a bivariate polynomial
function get_conv_rad(A::Generic.Poly{T}) where {T}
    
    # Return the minimum radius found for all coefficients series
    return minimum(map(get_conv_rad, coefficients(A)))
end

# Returns the approximate convergence radius by linear regression for a series
function get_conv_rad(A::Generic.AbsSeries{T}) where {T}
    # Prepare coefficients
    (powers, coeffs) = get_coeffs(A)

    if (length(coeffs) < 2)
        return 10000.0
    end

    # Generate linear regression model
    model = lm(@formula(y ~ x), DataFrame(x=powers, y=coeffs))

    # Get slope equal to ln(radius)
    radius = coef(model)[2]

    return exp(radius)
end

# In other cases answer inf
function get_conv_rad(_A::U) where {U}
    return 100000.0
end

# Returns two vectors (coeff exponent, log(|coeff|)) to be used for linear regression
function get_coeffs(A::Generic.AbsSeries{T}) where {T}
    coeffs = CFloatElem[]
    powers = CFloatElem[]
    for i in 1:(A.length - 1)
        c = coeff(A, i)
        if c != 0
            push!(coeffs, log(abs(c)))
            push!(powers, i)
        end
    end

    return (powers, coeffs)
end

# Returns a random Hensel setup with polynomials g and h of random degree <= m in y and <= l in x.
function random_hensel_setup(m::Int, l::Int)
    QX, _ = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
    QXY, _ = polynomial_ring(QX,"y")

    F = 0
    g = 0
    h = 0
    s = 0
    t = 0

    while true
        gi = random_polynomial(AbstractAlgebra.QQ, QX, QXY, m, l; monic=true)
        hi = random_polynomial(AbstractAlgebra.QQ, QX, QXY, m, l; monic=true)
        F = gi*hi
    
        g = biv_truncate(gi, 1)
        h = biv_truncate(hi, 1)
        s,t = cofactors(g, h, AbstractAlgebra.QQ, QX, QXY)

        check1 = s*g + t*h == 1
        check2 = biv_truncate(F - g*h, 1) == 0

        if check1 && check2 
            break
        end

    end

    return (F,g,h,s,t)
end


function generate_multiple_examples(y_degree::Int, x_degree::Int, number::Int)
    examples = []
    for _ in 1:number
        push!(examples, random_hensel_setup(y_degree, x_degree))
    end
    return examples
end


function get_distance(A::Generic.Poly{T}, B::Generic.Poly{U}, l::Int, rho::Float64) where {T,U}
    DX, _ = power_series_ring(DFloat, l, "x"; model=:capped_absolute)
    DXY, _ = polynomial_ring(DX,"y")

    Ad = to_other_poly(A, DFloat, DX, DXY)
    Bd = to_other_poly(B, DFloat, DX, DXY)

    return biv_norm(Ad - Bd, rho, 1.0)
end