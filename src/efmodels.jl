
#################################################
#
#   Type Hierarchy
#
#################################################

# exponential family models
abstract EFModel{VForms}

# exponential family distribution, i.e. p(x; θ)
abstract EFDistribution{F<:VariateForm} <: EFModel{Tuple{F}}

typealias EFUnivariateDistribution EFDistribution{Univariate}
typealias EFMultivariateDistribution EFDistribution{Multivariate}

# exponential family conditional distribution, i.e. p(x | y; θ)
abstract EFCondDistribution{Fx<:VariateForm, Fy<:VariateForm} <: EFModel{Tuple{Fx,Fy}}


#################################################
#
#   Generic/fallback functions
#
#################################################

## for EFModel

@req_method nsamples(em::EFModel, x)

@req_method constbdf(em::EFModel)
@req_method logpartition(em::EFModel)

@req_method inner!(r::AbstractArray, em::EFModel, x)

inner(em::EFModel, θ, x) = inner!(Array(Float64, nsamples(em, x)), em, x)

function logupdf!(r::AbstractArray, em::EFModel, x)
    if constbdf(em)
        b = logbdf(em)
        inner!(r, em, x)
        for i = 1:length(r)
            @inbounds r[i] += b
        end
    else
        inner!(r, em, x)
        for i = 1:length(r)
            @inbounds r[i] += logbdf(em, x)
        end
    end
    r
end
logupdf(em::EFModel, x) = logupdf!(Array(Float64, nsamples(em, x)), em, x)

function logpdf!(r::AbstractArray, em::EFModel, x)
    a = logpartition(em)
    logupdf!(r, em, x)
    for i = 1:length(r)
        @inbounds r[i] -= a
    end
    r
end
logpdf(em::EFModel, x) = logpdf!(Array(Float64, nsamples(em, x)), em, x)


# EFUnivariateDistribution

nsamples(d::EFUnivariateDistribution, x::Number) = 1
nsamples(d::EFUnivariateDistribution, x::AbstractArray) = length(x)

@req_method inner(d::EFUnivariateDistribution, x::Float64)
inner(d::EFUnivariateDistribution, x::Real) = inner(d, f64(x))

@req_method logbdf(d::EFUnivariateDistribution, x::Float64)
logbdf(d::EFUnivariateDistribution, x::Real) = logbdf(d, x)

logupdf(d::EFUnivariateDistribution, x::Float64) = inner(d, x) + logbdf(d, x)
logupdf(d::EFUnivariateDistribution, x::Real) = logupdf(d, f64(x))

logpdf(d::EFUnivariateDistribution, x::Float64) = logupdf(d, x) - logpartition(d)
logpdf(d::EFUnivariateDistribution, x::Real) = logpdf(d, f64(x))

function inner!(r::AbstractArray, d::EFUnivariateDistribution, x::AbstractArray)
    n = nsamples(d, x)
    length(r) == n || throw(DimensionMismatch("Inconsistent sizes of r and x."))
    for i = 1:length(x)
        @inbounds r[i] = inner(d, x[i])
    end
    r
end

# EFMultivariateDistribution

nsamples(d::EFMultivariateDistribution, x::AbstractVector) = 1
nsamples(d::EFMultivariateDistribution, x::AbstractMatrix) = size(x, 2)
