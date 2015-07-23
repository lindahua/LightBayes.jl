### Generic functions

### fallback generic functions

constbdf(d::Normal) =
    throw(MethodError(constbdf, (typeof(d),)))

logpartition(d::ExponentialFamily) =
    throw(MethodError(logpartition, (typeof(d),)))

inner!(r::AbstractArray, d::EFDistribution, x) =
    throw(MethodError(inner!, (typeof(r), typeof(d), typeof(x))))

inner(d::EFDistribution, x) = inner!(Array(Float64, nsamples(d, x)), d, x)

function logupdf!(r::AbstractArray, d::EFDistribution, x)
    if constbdf(d)
        b = logbdf(d)
        inner!(r, d, x)
        for i = 1:length(r)
            @inbounds r[i] += b
        end
    else
        inner!(r, d, x)
        for i = 1:length(r)
            @inbounds r[i] += logbdf(r, x)
        end
    end
    r
end
logupdf(d::EFDistribution, x) = logupdf!(Array(Float64, nsamples(d, x)), d, x)

function logpdf!(r::AbstractArray, d::EFDistribution, x)
    a = logpartition(d)
    logupdf!(r, d, x)
    for i = 1:length(r)
        @inbounds r[i] -= a
    end
    r
end
logpdf(d::EFDistribution, x) = logpdf!(Array(Float64, nsamples(d, x)), d, x)


# EFUnivariateDistribution

nsamples(d::EFUnivariateDistribution, x::Number) = 1
nsamples(d::EFUnivariateDistribution, x::AbstractArray) = length(x)

inner(d::EFUnivariateDistribution, x::Float64) =
    throw(MethodError(inner, (typeof(d), Float64)))

inner(d::EFUnivariateDistribution, x::Real) =
    inner(d, _fp(x))

logbdf(d::EFUnivariateDistribution, x::Float64) =
    throw(MethodError(logbdf, (typeof(d), x)))

logbdf(d::EFUnivariateDistribution, x::Real) =
    logbdf(d, x)

logupdf(d::EFUnivariateDistribution, x::Float64) =
    inner(d, x) + logbdf(d, x)

logpdf(d::EFUnivariateDistribution, x::Float64) =
    logupdf(d, x) - logpartition(d)

logupdf(d::EFUnivariateDistribution, x::Real) = logupdf(d, _fp(x))
logpdf(d::EFUnivariateDistribution, x::Real) = logpdf(d, _fp(x))

function inner!(r::AbstractArray, d::EFUnivariateDistribution, x::AbstractArray)
    n = length(x)
    length(r) == n || throw(DimensionMismatch("Inconsistent sizes of r and x."))
    for i = 1:length(x)
        @inbounds r[i] = inner(d, x[i])
    end
    r
end

# EFMultivariateDistribution

nsamples(d::EFMultivariateDistribution, x::AbstractVector) = 1
nsamples(d::EFMultivariateDistribution, x::AbstractMatrix) = size(x, 2)

# TODO: more methods for multivariates


### Specific exponential family distributions

module EFD

using Distributions
import Base: convert, show
import Distributions: nsamples, logpdf, log2π

import BayesBase:
    EFUnivariateDistribution,
    inner, logpartition, logbdf, constbdf,
    inner!, logupdf!, logupdf,
    _fp

include(joinpath("efd", "normal.jl"))

end
