immutable EFNormal <: EFUnivariateDistribution
    h::Float64      # = μ / σ^2
    J::Float64      # = 1 / σ^2

    EFNormal(h::Real) = new(f64(h), 1.0)

    function EFNormal(h::Real, J::Real)
        J > 0 || error("J must be positive.")
        new(f64(h), f64(J))
    end
end

show(io::IO, d::EFNormal) = print(io, "EFNormal(h=$(d.h), J=$(d.J))")

# conversion

function convert(::Type{Normal}, d::EFNormal)
    σ2 = 1.0 / d.J
    Normal(σ2 * d.h, sqrt(σ2))
end

function convert(::Type{EFNormal}, d::Normal)
    J = 1.0 / abs2(d.σ)
    EFNormal(d.μ * J, J)
end

convert(::Type{Distribution}, d::EFNormal) = convert(Normal, d)
convert(::Type{EFDistribution}, d::Normal) = convert(EFNormal, d)


# interface functions

logpartition(d::EFNormal) = 0.5 * abs2(d.h) / d.J

constbdf(d::EFNormal) = true
logbdf(d::EFNormal) = 0.5 * (log(d.J) - log2π)
logbdf(d::EFNormal, x::Float64) = logbdf(d)

inner(d::EFNormal, x::Float64) = d.h * x - (d.J * abs2(x)) * 0.5
