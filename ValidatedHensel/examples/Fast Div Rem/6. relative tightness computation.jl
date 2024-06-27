using AbstractAlgebra
using Nemo
using ValidatedHensel


function full_tightness_test(y_deg::Int, min_x_deg::Int, max_x_deg::Int, rho::Float64, tau::Float64, iterations::Int)
    for x_deg in min_x_deg:5:max_x_deg
        println("Test for degree: ", x_deg)
        (Q_tight, R_tight, exact_Q, exact_R, failures) = compute_relative_tightness(y_deg, x_deg, rho, tau, iterations)
        open("results/tightness_comparison.txt", "a+") do file
            write(file, string(x_deg)*" "*string(Q_tight)*" "*string(R_tight)*" "*string(exact_Q)*" "*string(exact_R)*" "*string(failures)*"\n")
        end
    end
end

function compute_relative_tightness(y_deg::Int, x_deg::Int, rho::Float64, tau::Float64, iterations::Int)
    total_Q = 0
    exact_Q = 0
    total_R = 0
    exact_R = 0
    failures = 0

    for _ in 0:(iterations-1)
        try
            relQ, relR, exQ, exR = iter_relative_tightness(y_deg, x_deg, rho, tau)
            if exQ
                exact_Q += 1
            else
                total_Q += relQ
            end

            if exR
                exact_R += 1
            else
                total_R += relR
            end
        catch
            failures+=1
        end
    end


    total_Q /= (iterations - failures - exact_Q)
    total_R /= (iterations - failures - exact_R)
    return (Float64(total_Q, RoundUp), Float64(total_R, RoundUp), exact_Q, exact_R, failures)
end

function iter_relative_tightness(y_deg::Int, x_deg::Int, rho::Float64, tau::Float64)
    # Q[[X]][Y]/(x^1000) for exact computation
    QX, _ = power_series_ring(AbstractAlgebra.QQ, x_deg, "x"; model=:capped_absolute)
    QXY, _ = polynomial_ring(QX,"y")

    # R[[X]][Y]/(x^1000) with 64-bit of precision for approximate computation
    RX, _ = power_series_ring(RDF, x_deg, "x"; model=:capped_absolute)
    RXY, _ = polynomial_ring(RX,"y")

    Fe = random_polynomial(AbstractAlgebra.QQ, QX, QXY, y_deg, x_deg)
    Ge = random_polynomial(AbstractAlgebra.QQ, QX, QXY, Int(ceil(y_deg/2)), x_deg; monic=true)


    # Exact computation
    (Qe,Re) = fast_div_rem(Fe, Ge, x_deg)

    # Sanity check
    if Fe - Ge*Qe - Re != 0
        error("Error when checking correctness of computation for \nF = " * string(Fe) *"\nG= " *string(Ge))
    end


    #Approximate computation
    Ff = to_other_poly(Fe, RDF, RX, RXY)
    Gf = to_other_poly(Ge, RDF, RX, RXY)

    (Qf, Rf, r_Q, r_R) = val_fast_div_rem(Ff, Gf, x_deg, rho, tau)

    
    #Comparison
    RRDF = typeof(ArbField(128)(1))
    RRX, _ = power_series_ring(ArbField(128), 10, "x"; model=:capped_absolute)
    RRY, _ = polynomial_ring(RRX, "y")

    Qef = to_other_poly(to_other_poly(Qe, RDF, RX, RXY), RRDF, RRX, RRY)
    Qff = to_other_poly(Qf, RRDF, RRX, RRY)
    q_norm= biv_norm(Qef-Qff, rho, tau)
    relative_Q = q_norm/r_Q
    
    Ref = to_other_poly(to_other_poly(Re, RDF, RX, RXY), RRDF, RRX, RRY)
    Rff = to_other_poly(Rf, RRDF, RRX, RRY)
    r_norm = biv_norm(Ref-Rff, rho, tau)
    relative_R = r_norm/r_R

    return (relative_Q, relative_R, r_Q == q_norm, r_R == r_norm)
end

full_tightness_test(20, 50, 100, 1.0, 1.0, 100)