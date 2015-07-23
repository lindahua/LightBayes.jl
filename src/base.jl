
abstract ExponentialFamily

abstract EFDistribution{F<:VariateForm} <: ExponentialFamily
abstract EFCondDistribution{Fx<:VariateForm, Fy<:VariateForm} <: ExponentialFamily

typealias EFUnivariateDistribution EFDistribution{Univariate}
typealias EFMultivariateDistribution EFDistribution{Multivariate}


### auxiliary functions

_fp(x::Real) = convert(Float64, x)


### fallback generic functions

logpartition(d::ExponentialFamily) =
    throw(MethodError(logpartition, (typeof(d),)))

# EFUnivariateDistribution

inner(d::EFUnivariateDistribution, x::Float64) =
    throw(MethodError(inner, (typeof(d), Float64)))

inner(d::EFUnivariateDistribution, x::Real) =
    inner(d, _fp(x))

logbdf(d::EFUnivariateDistribution, x::Float64) =
    throw(MethodError(logbdf, (typeof(d), x)))

logbdf(d::EFUnivariateDistribution, x::Real) =
    logbdf(d, x)

logupdf(d::EFUnivariateDistribution, x::Float64) =
    inner(d, x) + logbdf(x)

logpdf(d::EFUnivariateDistribution, x::Float64) =
    inner(d, x) - logpartition(d) + logbdf(x)

logupdf(d::EFUnivariateDistribution, x::Real) = logupdf(d, _fp(x))
logpdf(d::EFUnivariateDistribution, x::Real) = logpdf(d, _fp(x))
