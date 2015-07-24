
#################################################
#
#   Type Hierarchy
#
#################################################

# exponential family models
abstract EFModel{VForms}

# exponential family distribution, i.e. p(x; θ)
abstract EFDistribution{F<:VariateForm,S<:ValueSupport} <: EFModel{Tuple{F}}

typealias EFUnivariateDistribution{S<:ValueSupport} EFDistribution{Univariate,S}
typealias EFMultivariateDistribution{S<:ValueSupport} EFDistribution{Multivariate,S}


#################################################
#
#   Generic/fallback functions
#
#################################################

### for EFModel

@req_method nsamples(em::EFModel, x)

@req_method constbdf(em::EFModel)
@req_method logpartition(em::EFModel)

@req_method inner!(r::AbstractArray, em::EFModel, x)

inner(em::EFModel, θ, x) = inner!(Array(Float64, nsamples(em, x)), em, x)

function logupdf!(r::AbstractArray, em::EFModel, x)
    if constbdf(em)
        inner!(r, em, x)
        b = logbdf(em)
        if b != zero(b)
            for i = 1:length(r)
                @inbounds r[i] += b
            end
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


### for EFUnivariateDistribution

# nsamples

nsamples(d::EFUnivariateDistribution, x::Number) = 1
nsamples(d::EFUnivariateDistribution, x::AbstractArray) = length(x)

# inner

@req_method inner(d::EFUnivariateDistribution{Continuous}, x::Float64)
@req_method inner(d::EFUnivariateDistribution{Discrete}, x::Int)

inner(d::EFUnivariateDistribution{Continuous}, x::Real) = inner(d, Float64(x))
inner(d::EFUnivariateDistribution{Discrete}, x::Real) = inner(d, Int(x))

# inner!

function inner!(r::AbstractArray, d::EFUnivariateDistribution, x::AbstractArray)
    n = nsamples(d, x)
    length(r) == n || throw(DimensionMismatch("Inconsistent sizes of r and x."))
    for i = 1:length(x)
        @inbounds r[i] = inner(d, x[i])
    end
    r
end

# logbdf

@req_method logbdf(d::EFUnivariateDistribution{Continuous}, x::Float64)
@req_method logbdf(d::EFUnivariateDistribution{Discrete}, x::Int)

logbdf(d::EFUnivariateDistribution{Continuous}, x::Real) = logbdf(d, Float64(x))
logbdf(d::EFUnivariateDistribution{Discrete}, x::Real) = logbdf(d, Int(x))

# logupdf

logupdf(d::EFUnivariateDistribution{Continuous}, x::Float64) = inner(d, x) + logbdf(d, x)
logupdf(d::EFUnivariateDistribution{Discrete}, x::Int) = inner(d, x) + logbdf(d, x)

logupdf(d::EFUnivariateDistribution{Continuous}, x::Real) = logupdf(d, Float64(x))
logupdf(d::EFUnivariateDistribution{Discrete}, x::Real) = logupdf(d, Int(x))

# logpdf

logpdf(d::EFUnivariateDistribution{Continuous}, x::Float64) = logupdf(d, x) - logpartition(d)
logpdf(d::EFUnivariateDistribution{Discrete}, x::Int) = logupdf(d, x) - logpartition(d)

logpdf(d::EFUnivariateDistribution{Continuous}, x::Real) = logpdf(d, Float64(x))
logpdf(d::EFUnivariateDistribution{Discrete}, x::Real) = logpdf(d, Int(x))


### for EFMultivariateDistribution

nsamples(d::EFMultivariateDistribution, x::AbstractVector) = 1
nsamples(d::EFMultivariateDistribution, x::AbstractMatrix) = size(x, 2)
