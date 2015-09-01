using Distributions
using LightBayes
using Base.Test

# construct prior
d = 2
β0 = zeros(d)
σ0 = 10.0
κ0 = inv(σ0^2)
p0 = MvNormalCanon(β0, κ0)

@test length(p0) == d
@test mean(p0) == zeros(d)
@test var(p0) == fill(σ0^2, d)

# derive posterior
σ = 2.0
g = IsoGaussModel(σ)
X = randn(d, 10)

pp = posterior(p0, suffstats(g, X))

@test_approx_eq pp.h β0 + inv(σ^2) * vec(sum(X, 2))
@test_approx_eq pp.J.value κ0 + inv(σ^2) * size(X,2)
