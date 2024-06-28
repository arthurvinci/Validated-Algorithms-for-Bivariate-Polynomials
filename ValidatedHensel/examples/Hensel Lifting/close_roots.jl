using AbstractAlgebra
using Nemo
using ValidatedHensel

include("utils.jl")

################################################################################
#
#   Example where g and h are coprime because of rounding errors
#                 g = (y-epsilon)^m
#                 h = (y+epsilon)^m
#   
#                 F = g*h + y^(2m)*x^2
#
################################################################################

function test_close_roots()
    epsilon = 1.0
    for _ in 0:4
        epsilon /= 10
        for m in 1:50
            l = 1
            for _ in 0:9
                l *= 2
                try
                    (rv, error_Gv, error_Hv, error_G2, error_H2, rho)  = test_m_eps(m, epsilon, l)

                    open("results/close_roots.txt", "a+") do file
                        write(file, string(epsilon)*" "*string(m)*" "*string(l)*" "*string(rv)*" "*string(error_Gv)*" "*string(error_Hv)*" "*string(error_G2)*" "*string(error_H2)*" "*string(rho)*"\n")
                    end 
                catch
                    break;
                end
            end
        end

    end
end

function test_m_eps(m::Int, epsilon::Float64, l::Int)
    QX, x = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
    QXY, y = polynomial_ring(QX,"y")

    epsilon = AbstractAlgebra.QQ(epsilon)
    # Initialize the test
    g = (y - epsilon)^m
    h = (y + epsilon)^m
    F = g*h + y^(2m)* x^2
    (s,t) = cofactors(g, h, AbstractAlgebra.QQ, QX, QXY)

    # Compute the exact solution
    (ge, he, _, _) = fast_hensel_lifting(F,g,h,s,t,l)

    # Get the floating point equivalent of the setup
    RX, _ = power_series_ring(RDF, l, "x"; model=:capped_absolute)
    RXY, _ = polynomial_ring(RX,"y")

    Ff = to_other_poly(F, RDF, RX, RXY)
    gf = to_other_poly(g, RDF, RX, RXY)
    hf = to_other_poly(h, RDF, RX, RXY)
    sf = to_other_poly(s, RDF, RX, RXY)
    tf = to_other_poly(t, RDF, RX, RXY)

    # Compute the approximated maximum convergence radii
    (_, _, _, rho) = val_find_rho_log(Ff, gf, hf, sf, tf, l)    

    # Compute the validation with half of the radius
    (gv, hv, rv) = val_hensel_lifting(Ff, gf, hf, sf, tf, l, rho, 1.0)

    # Compute other algorithm output
    (g2, h2, _, _) = fast_hensel_lifting(Ff, gf, hf, sf, tf, l)

    # Compute error bounds with 128 bits precision to be rigorous
    RRDF = typeof(ArbField(128)(1))
    RRX, _ = power_series_ring(ArbField(128), l, "x"; model=:capped_absolute)
    RRY, _ = polynomial_ring(RRX, "y")

    Gef = to_other_poly(to_other_poly(ge, RDF, RX, RXY), RRDF, RRX, RRY)
    Hef = to_other_poly(to_other_poly(he, RDF, RX, RXY), RRDF, RRX, RRY)
    Gvf = to_other_poly(to_other_poly(gv, RDF, RX, RXY), RRDF, RRX, RRY)
    Hvf = to_other_poly(to_other_poly(hv, RDF, RX, RXY), RRDF, RRX, RRY)
    G2f = to_other_poly(to_other_poly(g2, RDF, RX, RXY), RRDF, RRX, RRY)
    H2f = to_other_poly(to_other_poly(h2, RDF, RX, RXY), RRDF, RRX, RRY)


    error_Gv = mag(biv_norm(Gef - Gvf, rho, 1.0))
    error_Hv = mag(biv_norm(Hef - Hvf, rho, 1.0))
    error_G2 = mag(biv_norm(Gef - G2f, rho, 1.0))
    error_H2 = mag(biv_norm(Hef - H2f, rho, 1.0))

    return (mag(rv), error_Gv, error_Hv, error_G2, error_H2, rho) 
end 

test_close_roots()