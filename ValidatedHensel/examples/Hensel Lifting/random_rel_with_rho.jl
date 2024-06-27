using AbstractAlgebra
using Nemo
using ValidatedHensel

include("utils.jl")




function iter_rho_incr(y_deg::Int, x_deg::Int)
    (F,g,h,s,t) = random_hensel_setup(y_deg, x_deg)

    # Compute exact solution first
    (ge, he, se_, te) = fast_hensel_lifting(F,g,h,s,t,x_deg)

    # Get the floating point equivalent of the setup
    RX, _ = power_series_ring(RDF, x_deg, "x"; model=:capped_absolute)
    RXY, _ = polynomial_ring(RX,"y")

    Ff = to_other_poly(F, RDF, RX, RXY)
    gf = to_other_poly(g, RDF, RX, RXY)
    hf = to_other_poly(h, RDF, RX, RXY)
    sf = to_other_poly(s, RDF, RX, RXY)
    tf = to_other_poly(t, RDF, RX, RXY)

    # Compute validation
    (gv, hv , rv) = val_hensel_lifting(Ff, gf, hf, sf, tf, x_deg, 1.0, 1.0)

    # Compute other output
    (g2, h2, _, _) = fast_hensel_lifting(Ff, gf, hf, sf, tf, x_deg)

    # Compute error bounds with 128 bits precision to be rigorous
    RRDF = typeof(ArbField(128)(1))
    RRX, _ = power_series_ring(ArbField(128), 10, "x"; model=:capped_absolute)
    RRY, _ = polynomial_ring(RRX, "y")

    Gef = to_other_poly(to_other_poly(ge, RDF, RX, RXY), RRDF, RRX, RRY)
    Hef = to_other_poly(to_other_poly(he, RDF, RX, RXY), RRDF, RRX, RRY)
    Gvf = to_other_poly(to_other_poly(gv, RDF, RX, RXY), RRDF, RRX, RRY)
    Hvf = to_other_poly(to_other_poly(hv, RDF, RX, RXY), RRDF, RRX, RRY)
    G2f = to_other_poly(to_other_poly(g2, RDF, RX, RXY), RRDF, RRX, RRY)
    H2f = to_other_poly(to_other_poly(h2, RDF, RX, RXY), RRDF, RRX, RRY)

    val_G_tightness = biv_norm(Gef - Gvf, rho, tau)/rv
    val_H_tightness = biv_norm(Hef - Hvf, rho, tau)/rv
    
    rel_alg_G_tightness = biv_norm(G2f - Gvf, rho, tau)/biv_norm(Gef - Gvf, rho, tau)
    rel_alg_H_tightness = biv_norm(H2f - Hvf, rho, tau)/biv_norm(Hef - Hvf, rho, tau)

    if isnan(mag(rel_alg_G_tightness)) || isnan(mag(rel_alg_G_tightness))
        println(biv_norm(Gef - Gvf, rho, tau))
        println(biv_norm(Hef - Hvf, rho, tau))
    end

    return (val_G_tightness, val_H_tightness, rel_alg_G_tightness, rel_alg_H_tightness)
end