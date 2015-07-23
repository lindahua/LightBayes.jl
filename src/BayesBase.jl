module BayesBase

using Distributions

import Base: convert, show
import Distributions: nsamples, logpdf, logpdf!

export
    # sub-modules
    EFD,

    # types
    ExponentialFamily,
    EFDistribution,
    EFUnivariateDistribution,
    EFMultivariateDistribution,

    # generic functions
    constbdf,
    inner, inner!,
    logpartition,
    logupdf, logupdf!,
    logbdf


## source files

include("base.jl")
include("efd.jl")


end # module
