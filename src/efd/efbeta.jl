immutable EFBeta <: EFUnivariateDistribution{Continuous}
    αm1::Float64
    βm1::Float64

    EFBeta() = new(0.0, 0.0)

    function EFBeta(αm1::Float64, βm1::Float64)
        (αm1 > -1.0 && βm1 > -1.0) ||
            error("Both αm1 and βm1 must be greater than -1.0.")
        new(αm1, βm1)
    end
    EFBeta(αm1::Real, βm1::Real) = EFBeta(Float64(αm1), Float64(βm1))
end

show(io::IO, d::EFBeta) = print(io, "EFBeta(αm1=$(d.αm1), βm1=$(d.βm1))")

# conversion

convert(::Type{Beta}, d::EFBeta) = Beta(d.αm1 + 1.0, d.βm1 + 1.0)
convert(::Type{EFBeta}, d::Beta) = EFBeta(d.α - 1.0, d.β - 1.0)

convert(::Type{Distribution}, d::EFBeta) = convert(Beta, d)
convert(::Type{EFDistribution}, d::Beta) = convert(EFBeta, d)

# interface functions

function logpartition(d::EFBeta)
    α = d.αm1 + 1.0
    β = d.βm1 + 1.0
    lgamma(α) + lgamma(β) - lgamma(α + β)
end

constbdf(d::EFBeta) = true
logbdf(d::EFBeta) = 0.0
logbdf(d::EFBeta, x::Float64) = 0.0

inner(d::EFBeta, x::Float64) = d.αm1 * log(x) + d.βm1 * log1p(-x)
