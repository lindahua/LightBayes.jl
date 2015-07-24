module EFD

using Distributions
import Base: convert, show
import Distributions: nsamples, logpdf, log2π

import BayesModels:
    EFUnivariateDistribution,
    inner, logpartition, logbdf, constbdf,
    inner!, logupdf!, logupdf,
    f64

include(joinpath("efd", "normal.jl"))

end
