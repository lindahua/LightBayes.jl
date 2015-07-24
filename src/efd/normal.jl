
immutable Normal <: EFUnivariateDistribution
    h::Float64      # =
    J::Float64

    Normal(h::Real) = new(f64(h), 1.0)

    function Normal(h::Real, J::Real)
        J > 0 || error("J must be positive.")
        new(f64(h), f64(J))
    end
end

show(io::IO, d::Normal) =
    print(io, "EFD.Normal(", d.h, ", ", d.J, ")")

# conversion

convert(::Type{Distribution}, d::Normal) =
    convert(Distributions.Normal, d)

function convert(::Type{Distributions.Normal}, d::Normal)
    σ2 = 1.0 / d.J
    Distributions.Normal(σ2 * d.h, sqrt(σ2))
end

function convert(::Type{Normal}, d::Distributions.Normal)
    J = 1.0 / abs2(d.σ)
    Normal(d.μ * J, J)
end

# interface functions

logpartition(d::Normal) = 0.5 * abs2(d.h) / d.J

constbdf(d::Normal) = true
logbdf(d::Normal) = 0.5 * (log(d.J) - log2π)
logbdf(d::Normal, x::Float64) = logbdf(d)

inner(d::Normal, x::Float64) = d.h * x - (d.J * abs2(x)) * 0.5
