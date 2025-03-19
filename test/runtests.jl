using LongMemory
using Test


@testset "ClassicEstimators" begin
    @test smean(ones(10,1)) == 1
    @test autocovariance(ones(10,1),1) ≈ [0.0]
    @test_throws ErrorException autocovariance(randn(10), 20)
    @test autocorrelation(collect(1:10),1) ≈ [1.0]
    @test autocorrelation(collect(1:10),1; flag=true) ≈ [1.0]
    @test_throws ErrorException autocorrelation(ones(10), 2)
    @test length(autocorrelation_plot(ones(10),1)) == 1
    @test length(log_variance_plot(randn(10,1); m = 4)) == 1
    @test length(log_variance_plot(randn(10,1); m = 4, slope=true, slope2=true)) == 1
    @test_throws ErrorException log_variance_plot(randn(2,1); m = 4)
    @test log_variance_est(collect(1:10); m = 4) ≈ 0.94596612
    @test_throws ErrorException log_variance_est(collect(1:10); m = 20)
    @test rescaled_range_est(collect(1:50)) ≈ -0.959937456
    @test_throws ErrorException rescaled_range_est(collect(1:10))
    @test length(rescaled_range_plot(collect(1:150))) == 1
    @test_throws ErrorException rescaled_range_plot(collect(1:10))
    @test rescaled_range(collect(1:30))[3] ≈ 33.0
    @test_throws ErrorException rescaled_range(collect(1:10))
    @test sstd(ones(10)) == 0.0
    @test sstdk(ones(10), 2) ≈ 0.0
end

@testset "DataExamples" begin
    @test NileData()[1,1] == 622
    @test NHTempData()[1,1] ≈ -0.901
    @test length(NileDataPlot()) == 4
    @test length(NHTempDataPlot()) == 4
    @test length(LMPlot(ones(100))) == 4
end

@testset "Forecasters" begin
    @test fi_ar_coefs(0.4,1) ≈ [0.4]
    @test fi_ar_coefs(0.2,1) ≈ [0.2]
    @test fi_ar_coefs(0.0,10) ≈ zeros(10,1)
    @test gregory_coeffs(3) ≈ [1/2; 1/12; 1/24]
    @test hitransform(hwfilter(ones(10))) ≈ ones(10)
    @test hitransform((zeros(10))) ≈ zeros(10)
    @test hitransform(ones(3)) ≈ [1; 1-1/2; 1-1/2-1/12]
    @test fi_forecast(ones(10), 1, 0)[11] ≈ 0.0
    @test fi_forecast(ones(10), 1, 0, 1)[11] ≈ 0.0
    @test length(fi_forecast_plot(ones(10), 1, 0)) == 1
    @test length(fi_forecast_plot(ones(10), 1, 0, 1)) == 1
    @test csa_forecast(zeros(10), 2, 1.5, 1.5) == zeros(12,1)
    @test csa_forecast(zeros(10), 2, 1.5, 1.5, 1)[1:12,1] == zeros(12)
    @test length(csa_forecast_plot(zeros(10), 2, 1.5, 1.5)) == 1
    @test length(csa_forecast_plot(zeros(10), 2, 1.5, 1.5, 1)) == 1
    @test LongMemory.Forecasters.csa_ma_coefs(1.5, 1.5, 2) ≈ [1.0; sqrt(1/2)]
    @test har_forecast(ones(100), 10)[101] ≈ 1.0
    @test length(har_forecast_plot(ones(100), 10)) == 1
    @test length(hwp_gen(10)) == 10
end

@testset "GeneratingFunctions" begin
    @test fracdiff(ones(10,1),0) == ones(10,1)
    @test fracdiff(ones(10,1),0.0) ≈ ones(10,1)
    @test fracdiff(ones(3,1),0.5) ≈ [1; 0.5; 0.375]
    @test fracdiff(ones(10), 1) == zeros(9)
    @test csadiff(zeros(10,1),1.4,1.4) ≈ zeros(10,1)
    @test length(csa_gen(10, 1.5, 1.5)) == 10
    @test length(csa_gen(10, 2, 1.5, 1.5)) == 10
    @test arfima_gen(10,1,[],0,[];σ=0) ≈ zeros(10,1)
    @test length(arfima_gen(10,1,[1],0,[];σ=0)) == 10
    @test length(arfima_gen(10,1,[1],0,[1];σ=0)) == 10
    @test arfi_gen(10,0,[],0;σ=0) ≈ zeros(10,1)
    @test arfi_gen(10,1,[],0;σ=0) ≈ ones(10,1)
    @test length(arfi_gen(10,1,[1],0;σ=0)) == 10
    @test length(fi(10,0.4)) == 10
    @test fi_gen(10,0;σ=0) ≈ zeros(10,1)
    @test length(edm_gen(10,0.4)) == 210
    @test length(sds_gen(10,0.4)) == 210    
    @test fi_survival_probs(3, 0)[:] == [1.0, 0.0, 0.0]
    @test hwfilter(zeros(10)) ≈ zeros(10)
    @test hwfilter((ones(3))) ≈ [1; 1+1/2; 1+1/2+1/3]

end

@testset "LongMemory.jl" begin


    @test fi_var_vals(1,0.4) ≈ [2.070098325296286]
    @test fi_var_vals(1,0.2) ≈ [1.0986855396043997]
    @test csa_var_vals(1, 1.4, 1.4) ≈ [2.1130846015858644]
    @test fi_cor_vals(1, 0.4) ≈ [1]
    @test csa_cor_vals(1, 1.4, 1.4) ≈ [1]
    
    
end
