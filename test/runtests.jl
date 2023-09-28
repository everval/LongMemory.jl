using LongMemory
using Test

@testset "LongMemory.jl" begin
    @test fracdiff(ones(10,1),0) ≈ ones(10,1)
    # Write your tests here.
end
