# Gaussian models

##
# Gaussian distribution with isotropic covariance
# as a conjugate prior.
#
immutable IsoGaussPrior <: PriorModel
    β::Vector{Float64}
    κ::Float64
end

length(d::IsoGaussPrior) = length(d.β)

mean(d::IsoGaussPrior) = inv(d.κ) * d.β
var(d::IsoGaussPrior) = inv(d.κ)

# log-partition function
function logpar(d::IsoGaussPrior)
    q = length(d)
    invκ = inv(d.κ)
    invκ * hsqrnorm(d.β) + (q / 2) * log(2π * invκ)
end

##
# Gaussian likelihood model
#
#  x ~ N(θ, σ^2)
#
immutable IsoGaussModel <: LikelihoodModel
    σ::Float64    # the standard dev of observations
end

function rand!(g::IsoGaussModel, θ::AbstractVector, X::AbstractMatrix)
    m, n = size(X)
    σ = g.σ
    length(θ) == m || throw(DimensionMismatch())
    @inbounds for j = 1:n
        for i = 1:m
            X[i,j] = θ[i] + rand() * σ
        end
    end
    X
end

rand(g::IsoGaussModel, θ::AbstractVector, n::Integer) =
    rand!(g, θ, Array(Float64, length(θ), n))

function posterior(d::IsoGaussPrior,           # prior distribution
                   g::IsoGaussModel,           # likelihood model
                   X::AbstractMatrix,          # data set
                   sel::AbstractVector{Int})   # selected indices

    β = copy(d.β)
    κ = d.κ
    c = 1.0 / g.σ^2
    n = length(sel)
    for i in sel
        x = view(X, :, i)
        axpy!(c, x, β)
    end
    IsoGaussPrior(β, κ + n * c)
end

posterior(d::IsoGaussPrior, g::IsoGaussModel, X::AbstractMatrix) =
    posterior(d, g, X, 1:size(X,2))
