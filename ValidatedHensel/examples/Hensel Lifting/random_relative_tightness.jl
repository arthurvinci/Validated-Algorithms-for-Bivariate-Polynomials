using AbstractAlgebra
using Nemo
using ValidatedHensel

include("utils.jl")

function tightness_in_x_incr(y_deg::Int, min_x_deg::Int, max_x_deg::Int, iterations::Int)
    for x_deg in min_x_deg:5:max_x_deg
        println("Test for degree: ", x_deg)
        compute_relative_tightness(y_deg, x_deg, iterations, "results/x_incr_tightness_comparison.txt")
    end
end

function compute_relative_tightness(y_deg::Int, x_deg::Int, iterations::Int, filename::String)
    arr = String[]
    for i in 0:(iterations-1)
        println("Iteration: ", i)
        (rv, error_Gv, error_Hv, error_G2, error_H2, rho) = iter_relative_tightness(y_deg, x_deg)
        if rv != 0
            push!(arr, string(x_deg)*" "*string(rv)*" "*string(error_Gv)*" "*string(error_Hv)*" "*string(error_G2)*" "*string(error_H2)*" "*string(rho))
        end 
    end
    open(filename, "a+") do file
        for result in arr
            write(file, result*"\n")
        end
        
    end
end

function iter_relative_tightness(y_deg::Int, x_deg::Int)
    (F,g,h,s,t) = random_hensel_setup(y_deg, x_deg)

    # Compute exact solution first
    (ge, he, _, _) = fast_hensel_lifting(F,g,h,s,t,x_deg)

    # Get the floating point equivalent of the setup
    RX, _ = power_series_ring(RDF, x_deg, "x"; model=:capped_absolute)
    RXY, _ = polynomial_ring(RX,"y")

    Ff = to_other_poly(F, RDF, RX, RXY)
    gf = to_other_poly(g, RDF, RX, RXY)
    hf = to_other_poly(h, RDF, RX, RXY)
    sf = to_other_poly(s, RDF, RX, RXY)
    tf = to_other_poly(t, RDF, RX, RXY)


    try
        # Compute validation with correct rho
        (gv, hv , rv, rho) = val_find_rho_log(Ff, gf, hf, sf, tf, x_deg)

        # Compute other output
        (g2, h2, _, _) = fast_hensel_lifting(Ff, gf, hf, sf, tf, x_deg)

        # Compute error bounds with 128 bits precision to be rigorous
        RRDF = typeof(ArbField(128)(1))
        RRX, _ = power_series_ring(ArbField(128), x_deg, "x"; model=:capped_absolute)
        RRY, _ = polynomial_ring(RRX, "y")
        

        Gef = to_other_poly(to_other_poly(ge, RDF, RX, RXY), RRDF, RRX, RRY)
        Hef = to_other_poly(to_other_poly(he, RDF, RX, RXY), RRDF, RRX, RRY)
        Gvf = to_other_poly(gv, RRDF, RRX, RRY)
        Hvf = to_other_poly(hv, RRDF, RRX, RRY)
        G2f = to_other_poly(g2, RRDF, RRX, RRY)
        H2f = to_other_poly(h2, RRDF, RRX, RRY)

        error_Gv = mag(biv_norm(Gef - Gvf, rho, 1.0))
        error_Hv = mag(biv_norm(Hef - Hvf, rho, 1.0))
        error_G2 = mag(biv_norm(Gef - G2f, rho, 1.0))
        error_H2 = mag(biv_norm(Hef - H2f, rho, 1.0))

        return (mag(rv), error_Gv, error_Hv, error_G2, error_H2, rho) 
    catch
        println("FAILED")
        return (0, 0, 0, 0, 0, 0)
    end
    
end


function same_degree_examples(min_degrees::Int, max_degrees::Int, iterations::Int)
    for i in min_degrees:max_degrees 
        
        sum_r = 0
        sum_rho = 0
        sum_g_incr = 0
        sum_h_incr = 0
        sum_rel_gn = 0
        sum_rel_hn = 0
        sum_rel_gf = 0
        sum_rel_hf = 0
        
        for j in 1:iterations
            println("degree: ", i, " iteration: ", j)
            # Sometimes the test fails for very specific cases (not due to the validation) so we loop until it works
            while true
                try
                    (r, rho, g_incr, h_incr, gn_dist, hn_dist, gfast_dist, hfast_dist) = random_hensel_test(i,i)
                    sum_r += r
                    sum_rho += rho
                    sum_g_incr += g_incr
                    sum_h_incr += h_incr
                    sum_rel_gn += r/gn_dist
                    sum_rel_hn += r/hn_dist
                    sum_rel_gf += r/gfast_dist
                    sum_rel_hf += r/hfast_dist
                    break
                catch
                end
            end
        end
        mean_r = sum_r/iterations
        mean_rho = rho/iterations
        mean_g_incr = sum_g_incr/iterations
        mean_h_incr = sum_h_incr/iterations
        mean_rel_gn = sum_rel_gn/iterations
        mean_rel_hn = sum_rel_hn/iterations
        mean_rel_gf = sum_rel_gf/iterations
        mean_rel_hf = sum_rel_hf/iterations

        open("results/hensel_same_degree_examples.txt", "a") do file
            write(file, string(i)*" "*string(mean_r)*" "*string(mean_rho)*" "*string(mean_g_incr)*" "*string(mean_h_incr)*" "*string(mean_rel_gn)*" "*string(mean_rel_hn)*" "*string(mean_rel_gf)*" "*string(mean_rel_hf))
        end
    end
end

function random_hensel_test(y_degree::Int, x_degree::Int)

    (F,g,h,s,t) = random_hensel_setup(y_degree, x_degree)

    # Compute exact solution first
    (ge, he, _, _) = fast_hensel_lifting(F,g,h,s,t,x_degree)

    # Get the floating point equivalent of the setup
    RX, _ = power_series_ring(RDF, x_degree, "x"; model=:capped_absolute)
    RXY, _ = polynomial_ring(RX,"y")

    Ff = to_other_poly(F, RDF, RX, RXY)
    gf = to_other_poly(g, RDF, RX, RXY)
    hf = to_other_poly(h, RDF, RX, RXY)
    sf = to_other_poly(s, RDF, RX, RXY)
    tf = to_other_poly(t, RDF, RX, RXY)

    # Compute a numeric solution
    (gn, hn, r, g_incr, h_incr, rho) = get_numeric_value_and_rhos(Ff, gf, hf, sf, tf, x_degree)

    # Compute distances to exact solution
    gn_dist = get_distance(ge, gn, x_degree, rho)
    hn_dist = get_distance(he, hn, x_degree, rho)

    # Compare with fast_numeric unrigorous approach
    (gfast, hfast, _, _) = fast_hensel_lifting(Ff, gf, hf, sf, tf, x_degree)
    gfast_dist = get_distance(ge, gfast, x_degree, rho)
    hfast_dist = get_distance(he, hfast, x_degree, rho)

    return (r, rho, g_incr, h_incr, gn_dist, hn_dist, gfast_dist, hfast_dist)
end

tightness_in_x_incr(20, 10, 50, 100)