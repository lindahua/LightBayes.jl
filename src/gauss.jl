# Gaussian models

### For prior IsoNormalCanon

function logpar(d::IsoNormalCanon)
    q = length(d)
    h = d.h
    σ2 = inv(d.J.value)
    (q * log(2π * σ2) + σ2 * vecnorm(h)^2) * 0.5
end

immutable IsoNormalCanonSuffStats
    h::Vector{Float64}  # add to h
    κ::Float64          # add to κ
end

length(ss::IsoNormalCanonSuffStats) = length(ss.h)

function posterior(d::IsoNormalCanon, ss::IsoNormalCanonSuffStats)
    length(d) == length(ss) || throw(DimensionMismatch())
    MvNormalCanon(d.h + ss.h, d.J.value + ss.κ)
end


##
# Gaussian likelihood model
#
#  x ~ N(θ, σ^2)
#
immutable IsoGaussModel <: LikelihoodModel
    σ::Float64    # the standard dev of observations
end

withparams(g::IsoGaussModel, μ::Vector{Float64}) =
    MvNormal(μ, g.σ)

function suffstats(g::IsoGaussModel, X::AbstractMatrix, inds::AbstractVector{Int})
    d, n = size(X)
    h = zeros(d)
    c = 1.0 / g.σ^2
    for i in inds
        x = view(X, :, i)
        axpy!(c, x, h)
    end
    IsoNormalCanonSuffStats(h, n * c)
end

suffstats(g::IsoGaussModel, X::AbstractMatrix) = suffstats(g, X, 1:size(X,2))
