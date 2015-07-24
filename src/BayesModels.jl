module BayesModels

using Distributions

import Base: convert, show
import Distributions: nsamples, logpdf, logpdf!
import Distributions: log2π, logistic, logit, log1pexp

export
    # efmodels
    EFModel,
    EFDistribution,
    EFUnivariateDistribution,
    EFMultivariateDistribution,

    # efd
    EFBernoulli,    # Bernoulli distribution
    EFNormal,       # Normal distribution

    constbdf,               # test whether the base density is constant
    inner, inner!,          # compute the inner product
    logpartition,           # compute the log-partition function
    logupdf, logupdf!,      # compute the logarithm of unnormalized pdf
    logbdf                  # compute the base density function


## source files

include("common.jl")
include("efmodels.jl")
include("efd.jl")

end # module
