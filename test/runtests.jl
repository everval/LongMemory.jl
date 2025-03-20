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
    @test length(rescaled_range_plot(collect(1:150); k=100, slope = true, slope2 = true)) == 1
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
    @test csa_ma_coefs(1.5, 1.5, 2) ≈ [1.0; sqrt(1/2)]
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

@testset "LogPeriodEstimators"  begin
    @test length(periodogram_plot(ones(10); slope = true)) == 1
    @test_throws ErrorException gph_est(ones(10); m = 0.1, l = 0.5)
    @test gph_est_variance(ones(10); m = 0.1, br = 0) ≈ π^2/24
    @test gph_est_variance(ones(10); m = 0.1, br = 1) ≈ π^2/24 * (9/4)
    @test gph_est_variance(ones(10); m = 0.1, br = 2) ≈ π^2/24 * (3.52)
    @test gph_est_variance(ones(10); m = 0.1, br = 3) ≈ π^2/24 * (4.79)
    @test gph_est_variance(ones(10); m = 0.1, br = 4) ≈ π^2/24 * (6.06)
    @test gph_est_variance(ones(10); m = 0.1, br = 5) ≈ π^2/24 * (6.06)
    @test gph_est_variance(ones(10); m = 0.1, br = -1) ≈ π^2/24 
    @test gph_est_variance(10; m = 0.1, br = 0) ≈ π^2/24
    @test gph_est_variance(10; m = 0.1, br = 1) ≈ π^2/24 * (9/4)
    @test gph_est_variance(10; m = 0.1, br = 2) ≈ π^2/24 * (3.52)
    @test gph_est_variance(10; m = 0.1, br = 3) ≈ π^2/24 * (4.79)
    @test gph_est_variance(10; m = 0.1, br = 4) ≈ π^2/24 * (6.06)
    @test gph_est_variance(10; m = 0.1, br = 5) ≈ π^2/24 * (6.06)
    @test gph_est_variance(10; m = 0.1, br = -1) ≈ π^2/24 
    @test LongMemory.LogPeriodEstimators.whittle_llk(0.4,collect(1:10)) ≈ 1.8053316105447856
    @test_throws ErrorException LongMemory.LogPeriodEstimators.whittle_llk(0.4,ones(10); m = 0.1, l = 0.5)
    @test whittle_est(collect(1:10)) ≈ 0.7496966402176871
    @test whittle_est_variance(collect(1:10)) ≈ 0.041666666666666664
    @test whittle_est_variance(10) ≈ 0.041666666666666664
    @test_throws ErrorException LongMemory.LogPeriodEstimators.exact_whittle_llk(0.4,collect(1:10); m = 0.1, l = 0.5)
    @test exact_whittle_est(collect(1:10)) ≈ 1.8101235290375945
    @test exact_whittle_est_variance(collect(1:10)) == whittle_est_variance(collect(1:10))
    @test exact_whittle_est_variance(10) == whittle_est_variance(10)
end

@testset "ParametricEstimators" begin
    @test LongMemory.ParametricEstimators.my_toeplitz(collect(1:5))[1,:] == 1:5
    @test fi_var_vals(1,0.4) ≈ [2.070098325296286]
    @test fi_var_vals(1,0.2) ≈ [1.0986855396043997]
    @test length(fi_var_vals(10, 0.4)) == 10
    @test size(fi_var_matrix(2, 0.1)) == (2,2)
    @test fi_var_matrix(10, 0.1)[1,1] ≈ 1.019494788
    @test fi_cor_vals(1, 0.4) ≈ [1]
    @test fi_llk(0.2, ones(10,1)) ≈ 1.0520643738861724
    @test length(fi_mle_est(rand(10,1))) == 2
    @test size(csa_var_matrix(5, 1.4, 1.4)) == (5,5)
    @test csa_var_vals(5, 1.4, 1.4)[1] ≈ 2.1130846015858644
    @test csa_cor_vals(1, 1.4, 1.4) ≈ [1]
    @test csa_cor_vals(5, 1.4, 1.4)[1,1] ≈ 1
    @test csa_llk(1.5,1.5,ones(10)) ≈ -0.08527518976186688
    @test length(csa_mle_est(ones(10))) == 3
    @test length(har_est(ones(50))) == 2
end

@testset "StructuralChanges" begin
    @test lm_change_test(ones(100)) == 15
    @test LongMemory.StructuralChanges.starred_var(ones(10))[1] ≈ 1.0
    @test LongMemory.StructuralChanges.simple_ols(ones(10),zeros(10)) ≈ 0.0
end
