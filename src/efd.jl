# All specific exponential family distributions

module EFD

import BayesBase:
    EFUnivariateDistribution,
    inner, logpartition, logbdf

include(joinpath("efd", "normal.jl"))

end
