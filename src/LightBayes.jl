module LightBayes

using StatsBase
using ArrayViews

import Base: length, mean, var
import Base.LinAlg: axpy!

export
    PriorModel,
    LikelihoodModel,
    IsoGaussPrior,
    IsoGaussModel,

    log_par,
    posterior

## source files

include("common.jl")
include("gauss.jl")

end # module
