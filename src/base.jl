
abstract ExponentialFamily

abstract EFDistribution{F<:VariateForm} <: ExponentialFamily
abstract EFCondDistribution{Fx<:VariateForm, Fy<:VariateForm} <: ExponentialFamily

typealias EFUnivariateDistribution EFDistribution{Univariate}
typealias EFMultivariateDistribution EFDistribution{Multivariate}


### auxiliary functions

_fp(x::Real) = convert(Float64, x)
