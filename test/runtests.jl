using LongMemory
using Test

@testset "LongMemory.jl" begin
    @test fracdiff(ones(10,1),0) == ones(10,1)
    @test fracdiff(ones(10,1),0.0) ≈ ones(10,1)
    @test fracdiff(ones(3,1),0.5) ≈ [1; 0.5; 0.375]
    @test csadiff(zeros(10,1),1.4,1.4) ≈ zeros(10,1)
    @test arfimagen(10,1,[],0,[];σ=0) ≈ zeros(10,1)
    @test arfigen(10,0,[],0;σ=0) ≈ zeros(10,1)
    @test arfigen(10,1,[],0;σ=0) ≈ ones(10,1)
    @test figen(10,0;σ=0) ≈ zeros(10,1)

end
