module LongMemory

export fracdiff, csadiff, csagen, edmgen, fi, figen, arfigen, arfimagen
include("GeneratingFunctions.jl")

export gph_est, gph_est_variance, whittle_est, exact_whittle_est, whittle_est_variance, periodogram
include("LogPeriodEstimators.jl")

end
