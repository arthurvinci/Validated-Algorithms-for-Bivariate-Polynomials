using AbstractAlgebra
using Nemo
using ValidatedHensel

function test_hensel_lifting_symbolic(k::Int, l::Int, iterations::Int=100)
    for _ in 0:iterations
        single_symbolic_hensel_test(k, l)
    end
    return true
end

function test_intermediate_hensel_lifting_symbolic(k::Int, l::Int, iterations::Int=100)
    for _ in 0:iterations
        single_symbolic_intermediate_hensel_test(k, l)
    end
    return true
end

function test_fast_hensel_lifting_symbolic(k::Int, l::Int, iterations::Int=100)
    for _ in 0:iterations
        single_symbolic_fast_hensel_test(k, l)
    end
    return true
end

function single_symbolic_hensel_test(k::Int, l::Int)
    # Q[[X]][Y]/(x^l) for exact computation
    QX, _ = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
    QXY, _ = polynomial_ring(QX,"y")

    F = 0
    g = 0
    h = 0
    s = 0
    t = 0

    while true
        gi = random_polynomial(AbstractAlgebra.QQ, QX, QXY, k, l)
        hi = random_polynomial(AbstractAlgebra.QQ, QX, QXY, k, l; monic=true)
        F = gi*hi
    
        g = biv_truncate(gi, 1)
        h = biv_truncate(hi, 1)
        s,t = cofactors(g, h, AbstractAlgebra.QQ, QX, QXY)

        check1 = s*g + t*h == 1
        check2 = biv_truncate(F - g*h, 1) == 0

        if check1 && check2 
            break
        else
            println(gi)
            println(hi)
        end

    end

    (G, H, S, T) = hensel_lifting(F, g, h, s, t, l)

    lift_check = F - G*H == 0
    cofactors_check = S*G + T*H == 1

    if !(lift_check && cofactors_check)
        error("Hensel Lifting failed for \nF: "*string(F)*"\ng: "*string(g)*"\nh: "*string(h)*"\ns: "*string(s)*"\nt: "*string(t))
    end
end


function single_symbolic_intermediate_hensel_test(k::Int, l::Int)
    # Q[[X]][Y]/(x^l) for exact computation
    QX, _ = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
    QXY, _ = polynomial_ring(QX,"y")

    F = 0
    g = 0
    h = 0
    s = 0
    t = 0

    while true
        gi = random_polynomial(AbstractAlgebra.QQ, QX, QXY, k, l)
        hi = random_polynomial(AbstractAlgebra.QQ, QX, QXY, k, l; monic=true)
        F = gi*hi

        g = biv_truncate(gi, 1)
        h = biv_truncate(hi, 1)
        s,t = cofactors(g, h, AbstractAlgebra.QQ, QX, QXY)

        check1 = s*g + t*h == 1
        check2 = biv_truncate(F - g*h, 1) == 0

        if check1 && check2 
            break
        else
            println(gi)
            println(hi)
        end

    end

    (G, H, S, T) = intermediate_hensel_lifting(F, g, h, s, t, l)

    lift_check = F - G*H == 0
    cofactors_check = S*G + T*H == 1

    if !(lift_check && cofactors_check)
        error("Intermediate Hensel Lifting failed for \nF: "*string(F)*"\ng: "*string(g)*"\nh: "*string(h)*"\ns: "*string(s)*"\nt: "*string(t))
    end
end


function single_symbolic_fast_hensel_test(k::Int, l::Int)
    # Q[[X]][Y]/(x^l) for exact computation
    QX, _ = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
    QXY, _ = polynomial_ring(QX,"y")

    F = 0
    g = 0
    h = 0
    s = 0
    t = 0

    while true
        gi = random_polynomial(AbstractAlgebra.QQ, QX, QXY, k, l)
        hi = random_polynomial(AbstractAlgebra.QQ, QX, QXY, k, l; monic=true)
        F = gi*hi

        g = biv_truncate(gi, 1)
        h = biv_truncate(hi, 1)
        s,t = cofactors(g, h, AbstractAlgebra.QQ, QX, QXY)

        check1 = s*g + t*h == 1
        check2 = biv_truncate(F - g*h, 1) == 0

        if check1 && check2 
            break
        else
            println(gi)
            println(hi)
        end

    end

    (G, H, S, T) = fast_hensel_lifting(F, g, h, s, t, l)

    lift_check = F - G*H == 0
    cofactors_check = S*G + T*H == 1

    if !(lift_check && cofactors_check)
        error("Fast Hensel Lifting failed for \nF: "*string(F)*"\ng: "*string(g)*"\nh: "*string(h)*"\ns: "*string(s)*"\nt: "*string(t))
    end
end