using LongMemory
using Test

@testset "LongMemory.jl" begin
    @test fracdiff(ones(10,1),0) == ones(10,1)
    @test fracdiff(ones(10,1),0.0) ≈ ones(10,1)
    @test fracdiff(ones(3,1),0.5) ≈ [1; 0.5; 0.375]
    @test csadiff(zeros(10,1),1.4,1.4) ≈ zeros(10,1)
    @test arfima_gen(10,1,[],0,[];σ=0) ≈ zeros(10,1)
    @test arfi_gen(10,0,[],0;σ=0) ≈ zeros(10,1)
    @test arfi_gen(10,1,[],0;σ=0) ≈ ones(10,1)
    @test fi_gen(10,0;σ=0) ≈ zeros(10,1)
    @test fi_ar_coefs(0.4,1) ≈ [0.4]
    @test fi_ar_coefs(0.2,1) ≈ [0.2]
    @test fi_ar_coefs(0.0,10) ≈ zeros(10,1)
    @test fi_var_vals(1,0.4) ≈ [2.070098325296286]
    @test fi_var_vals(1,0.2) ≈ [1.0986855396043997]
    @test csa_var_vals(1, 1.4, 1.4) ≈ [2.1130846015858644]
    @test fi_cor_vals(1, 0.4) ≈ [1]
    @test csa_cor_vals(1, 1.4, 1.4) ≈ [1]
    @test NileData()[1,1] == 622
    @test NHTempData()[1,1] ≈ -0.901
    @test hitransform(hwfilter(ones(10))) ≈ ones(10)
    @test gregory_coeffs(3) ≈ [1/2; 1/12; 1/24]
    @test hitransform((zeros(10))) ≈ zeros(10)
    @test hwfilter(zeros(10)) ≈ zeros(10)
    @test hwfilter((ones(3))) ≈ [1; 1+1/2; 1+1/2+1/3]
    @test hitransform(ones(3)) ≈ [1; 1-1/2; 1-1/2-1/12]
end
