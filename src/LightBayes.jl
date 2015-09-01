module LightBayes

using Distributions
using ArrayViews

import Base: length, mean, var, rand, rand!
import Base.LinAlg: axpy!
import Distributions: suffstats

export
    LikelihoodModel,
    IsoGaussModel,

    withparams,
    logpar,
    posterior,
    posterior_mode,
    posterior_mode!,
    posterior_rand,
    posterior_rand!

## source files

include("common.jl")
include("gauss.jl")

end # module
