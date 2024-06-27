using ApproximateRoot
using Test

include("test_app_root.jl")
include("test_numeric_app_root.jl")

@testset "ApproximateRoot.jl" begin
    @test TestAppRootV1()
    @test TestAppRootV2()
    @test TestNumericAppRootV1()
    @test TestNumericAppRootV2()
end

