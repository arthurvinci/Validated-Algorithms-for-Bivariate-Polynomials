using ValidatedHensel
using Test
using AbstractAlgebra

include("utils.jl")
include("test_fast_div_rem.jl")
include("test_hensel_lifting.jl")

@testset "Fast Div Rem" begin
    #@test test_fast_div_rem_symbolic(50, 50)
    #@test test_val_fast_div_rem(50, 50)
end

@testset "Hen*sel Lifting" begin
    @test test_hensel_lifting_symbolic(20, 20, 10)
    @test test_intermediate_hensel_lifting_symbolic(20, 20, 10)
    @test test_fast_hensel_lifting_symbolic(20, 20, 10)
end