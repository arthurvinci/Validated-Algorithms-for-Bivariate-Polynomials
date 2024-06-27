using Nemo
using AbstractAlgebra

precision(BigFloat)
eps(BigFloat)

FB = Nemo.AbstractAlgebra.Floats{BigFloat}()

function TestNumericAppRootV1()
    A, x = polynomial_ring(FB, "x")
    L, y = polynomial_ring(A,"y")
    F= 1.3*y^4+2.567*x^3*y^3+208.123*x^3*y^2+x^6

    phi=approximate_root(F, 1, 1)
    if (F-phi != 0)
        return false
    end
    
    phi=approximate_root(F, 4, 1)
    phi4 = phi^4

    if (degree(F-phi4)>2)
        return false
    end

    phi=approximate_root(F, 2, 1)
    phi2 = phi^2
    
    if (degree(F-phi2)>2)
        return false
    end

    return true
end

function TestNumericAppRootV2()
    A, x = polynomial_ring(FB, "x")
    L, y = polynomial_ring(A,"y")
    F=1.3*y^4+2.567*x^3*y^3+208.123*x^3*y^2+x^6

    phi=approximate_root(F, 1, 2)
    if (F-phi != 0)
        return false
    end
    
    phi=approximate_root(F, 4, 2)
    phi4 = phi^4

    if (degree(F-phi4)>2)
        return false
    end

    phi=approximate_root(F, 2, 2)
    phi2 = phi^2
    
    if (degree(F-phi2)>2)
        return false
    end

    return true
end