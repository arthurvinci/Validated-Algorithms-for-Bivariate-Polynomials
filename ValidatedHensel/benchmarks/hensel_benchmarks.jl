using AbstractAlgebra
using ValidatedHensel
using CPUTime

function same_degree_benchmark(max_degree::Int, nb_per_benchmark::Int)
    for i in 10:max_degree
        fast_total = 0
        medium_total = 0
        slow_total = 0

        for _ in 1:nb_per_benchmark
            (F,g,h,s,t) = hensel_setup(i, i)
            fast_total += @CPUelapsed fast_hensel_lifting(F,g,h,s,t,i)
            medium_total += @CPUelapsed intermediate_hensel_lifting(F,g,h,s,t,i)
            slow_total += @CPUelapsed hensel_lifting(F,g,h,s,t,i)
        end
        open("results/hensel_same_degree_bench.txt", "a") do file
            write(file, string(i)*" "*string(fast_total/nb_per_benchmark)*" "*string(medium_total/nb_per_benchmark)*" "*string(slow_total/nb_per_benchmark)*"\n")
        end
    end
end

function x_degree_grow_benchmark(y_degree::Int, max_x_degree::Int, nb_per_benchmark::Int)
    for i in 10:max_x_degree
        fast_total = 0
        medium_total = 0
        slow_total = 0
        for _ in 0:nb_per_benchmark
            (F,g,h,s,t) = hensel_setup(y_degree, i)
            fast_total += @CPUelapsed fast_hensel_lifting(F,g,h,s,t,i)
            medium_total += @CPUelapsed intermediate_hensel_lifting(F,g,h,s,t,i)
            slow_total += @CPUelapsed hensel_lifting(F,g,h,s,t,i)
        end
        open("results/hensel_x_grow_bench.txt", "a+") do file
            write(file, string(i)*" "*string(fast_total/nb_per_benchmark)*" "*string(medium_total/nb_per_benchmark)*" "*string(slow_total/nb_per_benchmark)*"\n")
        end
    end
end

function y_degree_grow_benchmark(max_y_degree::Int, x_degree::Int, nb_per_benchmark::Int)
    for i in 10:max_y_degree
        fast_total = 0
        medium_total = 0
        slow_total = 0
        for _ in 0:nb_per_benchmark
            (F,g,h,s,t) = hensel_setup(i, x_degree)
            fast_total += @CPUelapsed fast_hensel_lifting(F,g,h,s,t,x_degree)
            medium_total += @CPUelapsed intermediate_hensel_lifting(F,g,h,s,t,x_degree)
            slow_total += @CPUelapsed hensel_lifting(F,g,h,s,t,x_degree)
        end
        open("results/hensel_y_grow_bench.txt", "a+") do file
            write(file, string(i)*" "*string(fast_total/nb_per_benchmark)*" "*string(medium_total/nb_per_benchmark)*" "*string(slow_total/nb_per_benchmark)*"\n")
        end
    end
end


function hensel_setup(m::Int, l::Int)
    QX, _ = power_series_ring(AbstractAlgebra.QQ, l, "x"; model=:capped_absolute)
    QXY, _ = polynomial_ring(QX,"y")

    RX, _ = power_series_ring(RDF, l, "x"; model=:capped_absolute)
    RXY, _ = polynomial_ring(RX,"y")

    F = 0
    g = 0
    h = 0
    s = 0
    t = 0

    while true
        gi = random_polynomial(AbstractAlgebra.QQ, QX, QXY, m, l)
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

    F = to_other_poly(F, RDF, RX, RXY)
    g = to_other_poly(g, RDF, RX, RXY)
    h = to_other_poly(h, RDF, RX, RXY)
    s = to_other_poly(s, RDF, RX, RXY)
    t = to_other_poly(t, RDF, RX, RXY)

    return (F,g,h,s,t)
end


same_degree_benchmark(100, 1000)
y_degree_grow_benchmark(100, 20, 1000)
x_degree_grow_benchmark(20, 100, 1000)

