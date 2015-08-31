module LightBayes

using StatsBase
using ArrayViews

import Base: length, mean, var, rand, rand!
import Base.LinAlg: axpy!

export
    PriorModel,
    LikelihoodModel,
    IsoGaussPrior,
    IsoGaussModel,

    logpar,
    posterior

## source files

include("common.jl")
include("gauss.jl")

end # module
