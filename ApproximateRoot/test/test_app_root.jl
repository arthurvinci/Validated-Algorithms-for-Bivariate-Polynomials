using AbstractAlgebra

function TestAppRootV1()
    A, x = power_series_ring(AbstractAlgebra.GF(211), 33, "x")
    L, y = polynomial_ring(A,"y")
    F=y^4+x^3*y^3+209*x^3*y^2+x^6

    phi=approximate_root(F, 1, 1)
    if (F-phi != 0)
        return false
    end
    
    phi=approximate_root(F, 4, 1)

    if (degree(F-phi^4)>2)
        return false
    end

    phi=approximate_root(F, 2, 1)
    
    if (degree(F-phi^2)>2)
        return false
    end

    return true
end

function TestAppRootV2()
    A, x = power_series_ring(AbstractAlgebra.GF(211), 33, "x")
    L, y = polynomial_ring(A,"y")
    F=y^4+x^3*y^3+209*x^3*y^2+x^6

    phi=approximate_root(F, 1, 2)
    if (F-phi != 0)
        return false
    end

    phi=approximate_root(F, 4, 2)

    if (degree(F-phi^4)>2)
        return false
    end

    phi=approximate_root(F, 2, 2)

    if (degree(F-phi^2)>2)
        return false
    end

    return true
end