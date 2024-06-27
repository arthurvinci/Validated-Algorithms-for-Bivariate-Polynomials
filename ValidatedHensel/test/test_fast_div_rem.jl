using AbstractAlgebra
using Nemo
using ValidatedHensel

function test_fast_div_rem_symbolic(k::Int, l::Int)
    for _ in 0:100
        single_symbolic_fdr_test(k, l)
    end
    return true
end

function single_symbolic_fdr_test(k::Int, l::Int)
    # Q[[X]][Y]/(x^l) for exact computation
    QX, _ = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
    QXY, _ = polynomial_ring(QX,"y")

    A = random_polynomial(AbstractAlgebra.QQ, QX, QXY, k, l)
    B = random_polynomial(AbstractAlgebra.QQ, QX, QXY, k, l; monic=true)

    (Q,R) = fast_div_rem(A, B, l)

    if A - B*Q - R != 0
        error("Fast div rem test failed for \nA: "*string(A)*"\nB: "*string(B))
    end
end

function test_val_fast_div_rem(k::Int, l::Int)
    for _ in 0:100
        single_val_fdr_test(k, l)
    end
    return true
end

function single_val_fdr_test(k::Int, l::Int)
    rho = 0.1
    tau = 1.0

    QX, _ = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
    QXY, _ = polynomial_ring(QX,"y")

    A = random_polynomial(AbstractAlgebra.QQ, QX, QXY, k, l)
    B = random_polynomial(AbstractAlgebra.QQ, QX, QXY, A.length - 1, l; monic=true)

    (Q, R) = fast_div_rem(A, B, l)


    RX, _ = power_series_ring(RDF, l, "x"; model=:capped_absolute)
    RXY, _ = polynomial_ring(RX,"y")

    AF = to_other_poly(A, RDF, RX, RXY)
    BF = to_other_poly(B, RDF, RX, RXY)



    try
        (QF, RF, r_Q, r_R) = val_fast_div_rem(AF, BF, l, rho, tau)
        Q_prime = to_other_poly(Q, RDF, RX, RXY)
        actual_Q_bound = biv_norm(Q_prime - QF, rho, tau)
        given_Q_bound = mag(r_Q)
    
        R_prime = to_other_poly(R, RDF, RX, RXY)
        actual_R_bound = biv_norm(R_prime - RF, rho, tau)
        given_R_bound = mag(r_R)
    
        if given_Q_bound < actual_Q_bound || given_R_bound < actual_R_bound
            error("Computed bounds are too tight!")
        end
    catch e
        error("Validated fast_div_rem failed with error: "*string(e)*"\nA: "*string(A)*"\nB: "*string(B))
    end


end

