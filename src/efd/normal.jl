
immutable Normal <: EFUnivariateDistribution
    h::Float64
    J::Float64
end

logbdf(d::Normal, x::Float64) = 0.5 * (log(d.J) - log2π)
inner(d::Normal, x::Float64) = d.h * x - 0.5 * (d.J * abs2(x))
logpartition(d::Normal) = 0.5 * abs2(d.h) / d.J
