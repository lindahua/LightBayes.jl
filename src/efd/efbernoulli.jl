immutable EFBernoulli <: EFUnivariateDistribution{Discrete}
    θ::Float64      # = log(p / (1 - p))

    EFBernoulli() = new(0.0)
    EFBernoulli(θ::Real) = new(f64(θ))
end

show(io::IO, d::Bernoulli) = print(io, "EFBernoulli(θ=$(d.θ))")

# conversion

convert(::Type{Bernoulli}, d::EFBernoulli) = Bernoulli(logistic(d.θ))
convert(::Type{EFBernoulli}, d::Bernoulli) = EFBernoulli(logit(d.p))

convert(::Type{Distribution}, d::EFBernoulli) = convert(Bernoulli, d)
convert(::Type{EFDistribution}, d::Bernoulli) = convert(EFBernoulli, d)

# interface functions

logpartition(d::EFBernoulli) = log1pexp(d.θ)

constbdf(d::EFBernoulli) = true
logbdf(d::EFBernoulli) = 0.0
logbdf(d::EFBernoulli, x::Bool) = 0.0
logbdf(d::EFBernoulli, x::Int) = 0.0
logbdf(d::EFBernoulli, x::Real) = 0.0

inner(d::EFBernoulli, x::Bool) = ifelse(x, d.θ, 0.0)
inner(d::EFBernoulli, x::Int) = inner(d, Bool(x))
inner(d::EFBernoulli, x::Real) = inner(d, Bool(x))

logupdf(d::EFBernoulli, x::Bool) = inner(d, x)
logupdf(d::EFBernoulli, x::Int) = inner(d, Bool(x))
logupdf(d::EFBernoulli, x::Real) = inner(d, Bool(x))
