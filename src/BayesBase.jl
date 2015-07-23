module BayesBase

using Distributions

import Distributions: log2π
import Distributions: logpdf

export
    # sub-modules
    EFD,

    # types
    ExponentialFamily,
    EFDistribution,
    EFUnivariateDistribution,
    EFMultivariateDistribution,

    # generic functions
    inner,
    logpartition,
    logbdf


## source files

include("base.jl")
include("efd.jl")


end # module
