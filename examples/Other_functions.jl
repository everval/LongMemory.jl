using Distributions

# To check speed of betas against recursive
coefs = (beta.(p .+ (0:T-1), q) ./ beta(p, q)) .^ (1 / 2)