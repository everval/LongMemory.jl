using Distributions

function fi_survival_probs(N::Int, d)

   p = zeros(N, 1)
   p[1, 1] = 1
   for ii = 1:N-1
       p[ii+1, 1] = p[ii, 1] * (ii - 1 + d) / (ii + 1 - d)
   end

   return p

end

p = fi_survival_probs(10, 0.5)

T = 10

u = rand(10, 1)

lo = searchsortedlast.(Ref(p[:,1]), u, rev = true)

[u lo p]


function edmgen(T::Int, d; t=0.5, μ=0, σ=1)

   t = Int(round(T * t)) # Pre-sample

   p = fi_survival_probs(T+t, d) # Survival probabilities
   
   u = rand(T+t, 1) 

   tiempos = searchsortedlast.(Ref(p[:,1]), u, rev = true) # Random times

   x = zeros(T+t, 1) # Pre-allocate
   err = rand(Normal(μ,σ),T+t) # Error term

   for ii = 1:(T+t)
      x[ii:min((ii+tiempos[ii]-1),T+t),1] .= err[ii]
   end

   return x

end

function edmgen_Test(T::Int, d; t=0.5)

   t = Int(round(T * t)) # Pre-sample

   p = fi_survival_probs(T+t, d) # Survival probabilities
   
   u = rand(T+t, 1) 

   tiempos = searchsortedlast.(Ref(p[:,1]), u, rev = true) # Random times

   x = zeros(T+t, 1) # Pre-allocate
   err = ones(T+t,1) # Error term

   for ii = 1:(T+t)
      uxu = zeros(T+t,1) # Pre-allocate
      uxu[ii:min((ii+tiempos[ii]-1),T+t),1] .= err[ii] # Error term surviving
      x = x + uxu 
   end

   return x, tiempos, p, t

end