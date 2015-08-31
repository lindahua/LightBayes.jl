# common facilities

### value type conversion

f64(x::Real) = Float64(x)

### Common type hierarchy

##
#  The base type for all conjugate prior
#  distributions.
#
abstract PriorModel

##
#  The base type for all likelihood model
#
abstract LikelihoodModel

### Auxiliary functions

# half of squared L2-norm
hsqrnorm(x::AbstractArray) = vecnorm(x)^2 / 2
