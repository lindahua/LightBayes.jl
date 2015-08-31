using LightBayes
using Base.Test

# construct prior
d = 2
β0 = zeros(d)
σ0 = 10.0
κ0 = inv(σ0^2)
p0 = IsoGaussPrior(β0, κ0)

@test length(p0) == d
@test mean(p0) == zeros(d)
@test var(p0) == σ0^2

# derive posterior
σ = 2.0
g = IsoGaussModel(σ)
X = [1.0 0.0; -1.0 0.0; 1.0 0.0].'
@assert size(X,1) == d

pp = posterior(p0, g, X)

@test pp.β == β0 + inv(σ^2) * vec(sum(X, 2))
@test pp.κ == κ0 + inv(σ^2) * size(X,2)
