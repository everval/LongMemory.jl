using LongMemory
using Test
using Random
using LinearAlgebra

# Without Plots, the plotting functions explain how to enable them
@testset "Plots extension not loaded" begin
    @test Base.get_extension(LongMemory, :LongMemoryPlotsExt) === nothing
    @test_throws "requires Plots.jl" LMPlot(randn(100))
    @test_throws "requires Plots.jl" NileDataPlot()
    @test_throws "requires Plots.jl" fi_forecast_plot(randn(100), 5, 0.3)
    @test_throws "requires Plots.jl" periodogram_plot(randn(100))
    @test length(NileData().NileMin) == length(NileData().Year) # everything else works without Plots
end

using Plots

@testset "Plots extension loaded" begin
    @test Base.get_extension(LongMemory, :LongMemoryPlotsExt) !== nothing
    @test length(LMPlot(randn(100))) == 4
    @test_throws MethodError LMPlot("not a series") # wrong arguments give the usual MethodError
end

@testset "ClassicEstimators" begin
    @test smean(ones(10,1)) == 1
    @test smean(1:10) == 5.5 # ranges and other AbstractArrays are accepted
    @test autocovariance(ones(10,1),1) ≈ [0.0]
    @test_throws ErrorException autocovariance(randn(10), 20)
    @test autocorrelation(collect(1:10),1) ≈ [1.0]
    @test_throws ErrorException autocorrelation(ones(10), 2)
    @test length(autocorrelation_plot(ones(10),1)) == 1
    @test length(log_variance_plot(randn(100,1); m = 40)) == 1
    @test length(log_variance_plot(randn(100,1); m = 40, slope=true, slope2=true)) == 1
    @test_throws ErrorException log_variance_plot(randn(2,1); m = 4)
    @test_throws ErrorException log_variance_plot(randn(10,1); m = 4) # fewer than 3 points in the regression
    @test_throws ErrorException log_variance_est(collect(1:10); m = 4)
    @test_throws ErrorException log_variance_est(collect(1:10); m = 20)
    @test rescaled_range_est(collect(1:50)) ≈ 0.512566592618326
    @test_throws ErrorException rescaled_range_est(collect(1:10))
    @test length(rescaled_range_plot(collect(1:150))) == 1
    @test length(rescaled_range_plot(collect(1:150); k=100, slope = true, slope2 = true)) == 1
    @test_throws ErrorException rescaled_range_plot(collect(1:10))
    @test rescaled_range(collect(1:30))[3] ≈ 3.0 # R/S of 1:7 around its own mean: R = 6, S = 2
    @test rescaled_range(1:30) ≈ rescaled_range(collect(1:30))
    @test_throws ErrorException rescaled_range(collect(1:10))
    @test sstd(ones(10)) == 0.0
    @test sstdk(ones(10), 2) ≈ 0.0
    @test sstdk(collect(1.0:10), 2) ≈ 8.0 # all 5 complete blocks are used
    @test autocovariance(collect(1.0:5), 2) ≈ [2.0; 0.8]

    # Asymptotic distributions of the short memory tests reproduce the tabulated quantiles
    @test LongMemory.ClassicEstimators.brownian_bridge_range_cdf(0.809) ≈ 0.025 atol = 1e-3
    @test LongMemory.ClassicEstimators.brownian_bridge_range_cdf(1.862) ≈ 0.975 atol = 1e-3
    @test LongMemory.ClassicEstimators.brownian_bridge_variance_cdf(0.1869) ≈ 0.95 atol = 1e-3
    @test LongMemory.ClassicEstimators.bartlett_lrv(collect(1.0:5), 1) ≈ 2.0 + 2 * (1/2) * 0.8
    let r = lo_rs_test(randn(500))
        @test length(r) == 2
        @test 0 <= r[2] <= 1
    end
    let r = vs_test(randn(500); q = 0)
        @test length(r) == 2
        @test 0 <= r[2] <= 1
    end
    @test length(dfa_plot(randn(200); slope = true, slope2 = true)) == 1
    @test length(dfa_fluctuation(randn(200))[1]) == length(dfa_fluctuation(randn(200))[2])
    @test_throws ErrorException dfa_est(randn(30))
end

@testset "DataExamples" begin
    @test NileData().Year[1] == 622
    @test NileData().NileMin[1] ≈ 1157.0
    @test NHTempData().NHTemp[1] ≈ -0.901
    @test length(NileDataPlot()) == 4
    @test length(NHTempDataPlot()) == 4
    @test length(LMPlot(ones(100))) == 4
end

@testset "Forecasters" begin
    @test fi_ar_coefs(0.4,1) ≈ [0.4]
    @test fi_ar_coefs(0.2,1) ≈ [0.2]
    @test fi_ar_coefs(0.0,10) ≈ zeros(10,1)
    @test fi_ma_coefs(0.4, 3) ≈ [1; 0.4; 0.4 * 1.4 / 2]
    @test fi_ma_coefs(0.0, 3) ≈ [1; 0; 0]
    @test gregory_coeffs(3) ≈ [1/2; 1/12; 1/24]
    @test hitransform(hwfilter(ones(10))) ≈ ones(10)
    @test hitransform((zeros(10))) ≈ zeros(10)
    @test hitransform(ones(3)) ≈ [1; 1-1/2; 1-1/2-1/12]
    @test hitransform([1, 1, 1]) ≈ [1; 1-1/2; 1-1/2-1/12]
    @test hitransform([2.0]) == [2.0]
    @test gregory_coeffs(4)[4] ≈ 19/720
    @test fi_forecast(ones(10), 1, 0)[11] ≈ 1.0 # the sample mean is added back
    @test fi_forecast(ones(10), 1, 0, 1)[11] ≈ 1.0
    @test fi_forecast([1.0; -1.0], 1, 0.4)[3] ≈ 0.4 * (-1.0) + 0.4 * 0.6 / 2 * 1.0
    # bands are ± 2 standard errors of the k-step forecast error
    let xfor = fi_forecast(randn(50), 3, 0.4, 2.0)
        @test xfor[51:53, 3] - xfor[51:53, 1] ≈ 2 * 2.0 .* sqrt.(cumsum(fi_ma_coefs(0.4, 3) .^ 2))
    end
    @test length(fi_forecast_plot(ones(10), 1, 0)) == 1
    @test length(fi_forecast_plot(ones(10), 1, 0, 1)) == 1
    @test csa_forecast(zeros(10), 2, 1.5, 1.5) == zeros(12)
    @test csa_forecast(zeros(10), 2, 1.5, 1.5, 1)[1:12,1] == zeros(12)
    @test length(csa_forecast_plot(zeros(10), 2, 1.5, 1.5)) == 1
    @test length(csa_forecast_plot(zeros(10), 2, 1.5, 1.5, 1)) == 1
    @test csa_ma_coefs(1.5, 1.5, 2) ≈ [1.0; sqrt(1/2)]
    # CSA one-step forecast is Σ_{j≥1} ψ_j ε_{T+1-j}; for x = [1, -1] the shocks are e1 = 1, e2 = -1 - ψ_1
    let ψ = csa_ma_coefs(1.5, 1.5, 3)
        @test csa_forecast([1.0; -1.0], 1, 1.5, 1.5)[3] ≈ ψ[2] * (-1 - ψ[2]) + ψ[3]
    end
    @test har_forecast(ones(100), 10)[101] ≈ 1.0
    let x = sin.(1:40) .+ 0.01 .* (1:40), β = har_est(sin.(1:40) .+ 0.01 .* (1:40))[1], T = 40
        # one-step HAR forecast uses the regressors built from x[1:T]; works for T < 2*max(m)
        @test har_forecast(x, 1)[T+1] ≈ β[1] + β[2]*x[T] + β[3]*sum(x[T-4:T])/5 + β[4]*sum(x[T-21:T])/22
    end
    let x = sin.(1:200) .+ 0.1 .* cos.((1:200) .^ 2), xfor = har_forecast(sin.(1:200) .+ 0.1 .* cos.((1:200) .^ 2), 5; flag = true)
        σ = har_est(x)[2]
        @test xfor[201, 3] - xfor[201, 1] ≈ 2 * σ # one-step band is ± 2σ
        @test all(diff(xfor[201:205, 3] - xfor[201:205, 1]) .>= 0) # bands widen with the horizon
    end
    let m = [22, 1, 5]
        har_est(sin.(1:60); m = m)
        @test m == [22, 1, 5]
    end
    @test length(har_forecast_plot(ones(100), 10)) == 1
    @test length(hwp_gen(10)) == 10

    # ARFIMA forecasting
    @test arfima_ar_coefs(0.3, [], [], 20) ≈ fi_ar_coefs(0.3, 20)
    @test arfima_ma_coefs(0.3, [], [], 20) ≈ fi_ma_coefs(0.3, 20)
    let πc = arfima_ar_coefs(0.3, [0.5, -0.2], [0.4], 100), ψ = arfima_ma_coefs(0.3, [0.5, -0.2], [0.4], 100)
        a = [1; -πc[1:99]]
        @test [sum(a[1:k] .* reverse(ψ[1:k])) for k in 1:100] ≈ [1; zeros(99)] atol = 1e-10 # (1 - π(L))ψ(L) = 1
    end
    let y = sin.(1:50)
        @test arfima_forecast(y, 1, 0.0, [0.6], [])[51] ≈ smean(y) + 0.6 * (y[50] - smean(y)) # AR(1)
        @test arfima_forecast(y, 5, 0.3, [], []) ≈ fi_forecast(y, 5, 0.3)
        @test size(arfima_forecast(y, 5, 0.3, [0.2], [0.1], 1.0)) == (55, 3)
        @test length(arfima_forecast_plot(y, 5, 0.3, [0.2], [0.1])) == 1
        @test length(arfima_forecast_plot(y, 5, 0.3, [0.2], [0.1], 1.0)) == 1
    end
end

@testset "GeneratingFunctions" begin
    @test fracdiff(ones(10,1),0) == ones(10)
    @test fracdiff(ones(10,1),0.0) ≈ ones(10,1)
    @test fracdiff(ones(3,1),0.5) ≈ [1; 0.5; 0.375]
    @test fracdiff(ones(10), 1) == zeros(9)
    @test csadiff(zeros(10,1),1.4,1.4) ≈ zeros(10,1)
    @test length(csa_gen(10, 1.5, 1.5)) == 10
    @test length(csa_gen(10, 2, 1.5, 1.5)) == 10
    @test fracdiff(ones(3), -1) ≈ [1; 2; 3] # integration via the generic filter
    @test fracdiff(ones(3), 1//2) ≈ [1; 0.5; 0.375]
    @test arfima_gen(10,1,[],0,[];σ=0) ≈ ones(10) # μ is the mean of the series, also without AR terms
    @test arfima_gen(10,2,[0.5],0,[];σ=0) ≈ 2 .* ones(10)
    @test length(arfima_gen(10,1,[1],0,[];σ=0)) == 10
    @test length(arfima_gen(10,1,[1],0,[1];σ=0)) == 10
    @test arfi_gen(10,0,[],0;σ=0) ≈ zeros(10,1)
    @test arfi_gen(10,1,[],0;σ=0) ≈ ones(10,1)
    @test length(arfi_gen(10,1,[1],0;σ=0)) == 10
    @test length(fi_gen(10, 1)) == 10 # integer d: integration through the generic filter
    @test fi_gen(10,0;σ=0) ≈ zeros(10,1)
    @test fi_gen(10,0.4;μ=3,σ=0) ≈ 3 .* ones(10)
    @test csa_gen(10,1.4,1.4;μ=3,σ=0) ≈ 3 .* ones(10)
    @test edm_gen(10,0.4;μ=3,σ=0) ≈ 3 .* ones(10)
    @test length(edm_gen(10,0.4)) == 10 # pre-sample is discarded
    @test length(sds_gen(10,0.4)) == 10
    @test fi_survival_probs(3, 0)[:] == [1.0, 0.0, 0.0]
    @test hwfilter(zeros(10)) ≈ zeros(10)
    @test hwfilter((ones(3))) ≈ [1; 1+1/2; 1+1/2+1/3]
    @test length(davies_harte(0.5 .^ (0:99))) == 100
    @test length(davies_harte([1.0])) == 1
    @test_throws ErrorException davies_harte([1.0; 2.0]) # not a valid covariance
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
    @test exact_whittle_est(collect(1:10) .+ 5; demean = true) ≈ exact_whittle_est(collect(1:10); demean = true)
    @test exact_whittle_est_variance(collect(1:10)) == whittle_est_variance(collect(1:10))
    @test exact_whittle_est_variance(10) == whittle_est_variance(10)
    # bandwidth m = round(T^0.8) uses Fourier frequencies j = 1, ..., m (periodogram indices 2:m+1)
    @test LongMemory.LogPeriodEstimators.freq_band(1000, 0.8, 0) == 2:252
    @test LongMemory.LogPeriodEstimators.freq_band(10, 0.99, 0) == 2:6
    @test gph_est(1:100) ≈ gph_est(collect(1:100))
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
    @test csa_var_vals(5, 1.4, 1.4)[1] ≈ (1.4 + 1.4 - 1) / (1.4 - 1) # γ(0) = B(p, q-1)/B(p, q)
    @test csa_var_vals(1, 1.5, 3.0)[1] ≈ sum(abs2, csa_ma_coefs(1.5, 3.0, 100_000)) atol = 1e-6
    @test csa_cor_vals(1, 1.4, 1.4) ≈ [1]
    @test csa_cor_vals(5, 1.4, 1.4)[1,1] ≈ 1
    @test csa_llk(1.5,1.5,ones(10)) ≈ -0.08527518976186688
    @test length(csa_mle_est(ones(10))) == 3
    @test length(har_est(ones(50))) == 2

    # Durbin-Levinson reproduces the Cholesky-based Gaussian likelihood terms
    let x = sin.(1:200) .+ cos.(sqrt.(1:200)), V = fi_var_matrix(200, 0.3), C = cholesky(Symmetric(fi_var_matrix(200, 0.3)))
        xhat, v = durbin_levinson(fi_var_vals(200, 0.3), x)
        @test sum(log.(v)) ≈ logdet(C)
        @test sum((x .- xhat) .^ 2 ./ v) ≈ dot(x, C \ x)
    end
    @test_throws ErrorException durbin_levinson(ones(3), ones(5))

    # Exact simulation
    @test length(fi_exact_gen(100, 0.3)) == 100
    @test fi_exact_gen(10, 0.3; μ = 2, σ = 0) ≈ 2 .* ones(10)
    @test_throws ErrorException fi_exact_gen(100, 0.5)
    @test length(csa_exact_gen(100, 1.5, 1.6)) == 100
    @test_throws ErrorException csa_exact_gen(100, 1.5, 1.0)

    # ARFIMA internals: partial autocorrelations map to stationary AR polynomials
    @test LongMemory.ParametricEstimators.pacf_to_coefs([0.5]) ≈ [0.5]
    @test LongMemory.ParametricEstimators.pacf_to_coefs([0.5, 0.2]) ≈ [0.5 - 0.2 * 0.5, 0.2]
    @test LongMemory.ParametricEstimators.arfima_spec([π/2], 0.0, [], []) ≈ [1.0]
end

@testset "StructuralChanges" begin
    @test lm_change_test(ones(100)) == 15
    @test LongMemory.StructuralChanges.starred_var(ones(10))[1] ≈ 1.0
    @test LongMemory.StructuralChanges.simple_ols(ones(10),zeros(10)) ≈ 0.0
    let r = lm_change_test(randn(100); flag = true, nsim = 9)
        @test length(r) == 3
        @test 0 < r[3] <= 1
    end
    let r = qu_test(randn(500))
        @test length(r) == 2
        @test r[2] == [1.022, 1.155, 1.426]
    end
    @test_throws ErrorException qu_test(randn(500); ε = 0.1)
end

# Functions return vectors, not T×1 matrices, also for T×1 matrix input. Only outputs with several columns
# (forecasts with confidence bands, covariance matrices) are matrices.
@testset "Return shapes" begin
    xm = randn(100, 1)
    @test fracdiff(xm, 0) isa AbstractVector
    @test fracdiff(xm, 1) isa AbstractVector
    @test fracdiff(xm, 0.3) isa AbstractVector
    @test csadiff(xm, 1.5, 1.5) isa AbstractVector
    @test hwfilter(xm) isa AbstractVector
    @test csa_gen(10, 5, 1.5, 1.5) isa AbstractVector
    @test csa_gen(10, 1.5, 1.5) isa AbstractVector
    @test edm_gen(10, 0.4) isa AbstractVector
    @test sds_gen(10, 0.4) isa AbstractVector
    @test fi_gen(10, 0.4) isa AbstractVector
    @test arfima_gen(10, 0, [0.2], 0.1, [0.3]) isa AbstractVector
    @test arfi_gen(10, 0, [0.2], 0.1) isa AbstractVector
    @test fi_survival_probs(5, 0.3) isa AbstractVector
    @test autocovariance(xm, 3) isa AbstractVector
    @test autocorrelation(xm, 3) isa AbstractVector
    @test rescaled_range(xm) isa AbstractVector
    @test periodogram(xm)[1] isa AbstractVector
    @test fi_forecast(xm, 2, 0.3) isa AbstractVector
    @test csa_forecast(xm, 2, 1.5, 1.5) isa AbstractVector
    @test arfima_forecast(xm, 2, 0.3, [], []) isa AbstractVector
    @test har_forecast(xm, 2) isa AbstractVector
    @test hitransform(xm) isa AbstractVector
    @test har_est(xm)[1] isa AbstractVector
    @test size(fi_forecast(xm, 2, 0.3, 1.0)) == (102, 3)
end

# Statistical checks with fixed seeds. The tolerances are about four standard errors, so they do not depend on the
# particular random stream.
@testset "Monte Carlo checks" begin
    Random.seed!(20261015)

    x = fi_exact_gen(4000, 0.3)
    @test abs(gph_est(x) - 0.3) < 4 * sqrt(gph_est_variance(4000))
    @test abs(whittle_est(x) - 0.3) < 4 * sqrt(whittle_est_variance(4000))
    @test abs(exact_whittle_est(x) - 0.3) < 4 * sqrt(exact_whittle_est_variance(4000))
    @test abs(feasible_exact_whittle_est(5 .+ x) - 0.3) < 4 * sqrt(exact_whittle_est_variance(4000))
    @test abs(fi_mle_est(x[1:1000])[1] - 0.3) < 0.1
    @test abs(fi_mle_est(x[1:1000])[2] - 1.0) < 0.1
    @test abs(dfa_est(x) - 0.3) < 0.15

    # Unknown mean with d > 1/2: the feasible exact local Whittle estimator uses x[1] instead of the sample mean
    y = 5 .+ fi_gen(2000, 1.3)
    @test abs(feasible_exact_whittle_est(y; m = 0.65) - 1.3) < 4 * sqrt(1 / (4 * round(2000^0.65)))

    # ARFIMA(1,d,0) Whittle estimation
    z = arfima_gen(4000, 1.0, [0.5], 0.2, []; σ = 2.0)
    d, AR, MA, σ = arfima_whittle_est(z, 1, 0)
    @test abs(d - 0.2) < 0.15
    @test abs(AR[1] - 0.5) < 0.15
    @test isempty(MA)
    @test abs(σ - 2.0) < 0.15

    # Davies-Harte: sample autocovariances of the exact simulation
    S = zeros(2)
    for ii = 1:4000
        w = fi_exact_gen(20, 0.3)
        S .+= [w[1]^2; w[1] * w[20]]
    end
    @test S[1] / 4000 ≈ fi_var_vals(20, 0.3)[1] atol = 0.1
    @test S[2] / 4000 ≈ fi_var_vals(20, 0.3)[20] atol = 0.1

    # CSA MLE is deterministic and recovers σ
    c = csa_exact_gen(1000, 1.5, 1.6; σ = 2.0)
    @test csa_mle_est(c) == csa_mle_est(c)
    @test abs(csa_mle_est(c)[3] - 2.0) < 0.25

    # Short memory tests reject under long memory
    @test lo_rs_test(fi_exact_gen(4000, 0.4))[2] < 0.05
    @test vs_test(fi_exact_gen(4000, 0.4))[2] < 0.05
end

# Every `julia> ...` line in the docstrings of the exported functions must run.
@testset "Docstring examples" begin
    sandbox = Module()
    Core.eval(sandbox, :(using LongMemory))
    for name in names(LongMemory)
        obj = getfield(LongMemory, name)
        obj isa Function || continue
        doc = string(Base.Docs.doc(obj))
        for line in split(doc, '\n')
            sline = strip(line)
            startswith(sline, "julia> ") || continue
            ex = Meta.parse(sline[8:end])
            @test (Core.eval(sandbox, ex); true)
        end
    end
end
