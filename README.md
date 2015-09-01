# LightBayes

This is a Julia package that implements several light-weight Bayesian models, primarily for supporting the research on Bayesian inference and estimation.

**Note:** This package is still in the *experimental* stage, and is primarily for support internal research purpose. Also, it has not been registered at the official registry [METADATA.jl](https://github.com/JuliaLang/METADATA.jl).

## Setup

As it has not been officially registered, you may not use `Pkg.add` to install the package. Instead, you should check it out directly as

```bash
# enter the directory that hosts the Julia repositories
cd ~/.julia/v0.4

# clone the package
git clone https://github.com/lindahua/LightBayes.jl.git LightBayes

# run the tests to make sure it works
cd LightBayes
julia test/runtests.jl
```

## API

This package introduces an abstract type:

```julia
# The base type for all likelihood models,
# which connect the parameters with observations
abstract LikelihoodModel
```

Here, we assume that both the parameter space and observation space are always vector spaces.

#### Methods for Prior

Prior distributions are simply using the distributions in the [Distributions](https://github.com/JuliaStats/Distributions.jl) package. However, for those distributions that may serve as a *prior* here, we introduce additional methods:

```julia
# Let pri be a prior distribution

# compute the unnormalized log-pdf for given parameters
# θ can be either a single parameter or an array of
# multiple parameters
logupdf(pri, θ)
logupdf!(r, pri, θ)

# compute the log-partition value
logpar(pri)

# compute the posterior distribution, given sufficient
# statistics collected from observations
posterior(pri, sstats)

# find the mode of the posterior distribution
# This is useful for MAP estimation
posterior_mode(pri, sstats)
posterior_mode!(r, pri, sstats)

# sample from the posterior distribution
posterior_rand(pri, sstats)
posterior_rand!(r, pri, sstats)
```

#### Methods for Likelihood Model

```julia
# Let md be a likelihood model

# return a distribution, given parameter
d = withparams(md, θ)

# get the number of samples in X, w.r.t. md
n = nsamples(md, X)

# compute sufficient statistics, given observed data
ss = suffstats(md, X)
ss = suffstats(md, X, inds)

# to get a posterior given data, one can write
post = posterior(pri, suffstats(md, X))
```

## Models

Currently, we implement the following likelihood models

```julia
# x ~ N(θ, σ^2), where σ^2 is fixed a priori
# conjugate prior type: IsoNormal
immutable IsoGaussModel <: LikelihoodModel
    σ::Float64
end

```
