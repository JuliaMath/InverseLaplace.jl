# InverseLaplace
### Numerical inverse Laplace transform

Linux, OSX: [![Build Status](https://travis-ci.org/jlapeyre/InverseLaplace.jl.svg)](https://travis-ci.org/jlapeyre/InverseLaplace.jl)
[![Coverage Status](https://coveralls.io/repos/github/jlapeyre/InverseLaplace.jl/badge.svg)](https://coveralls.io/github/jlapeyre/InverseLaplace.jl)
[![codecov.io](http://codecov.io/github/jlapeyre/InverseLaplace.jl/coverage.svg?branch=master)](http://codecov.io/github/jlapeyre/InverseLaplace.jl?branch=master)

### ilt

Example. Evaluate the inverse Laplace transform of `1/s^2` at
`t=0.5`.
```julia
    ilt(s -> 1/s^2, 0.5)
```

The third argument gives the number of terms used.
Use more terms for BigFloat input
```julia
    ilt(s -> 1/s^2, BigFloat(1//2),80)
```

Defaults are `32` for `Float64` and `64` for `BigFloat`.


### gwr

`ilt` evaluates the function at complex values. The function `gwr` employs an algorithm
that evaluates the function only for real values. `gwr` is called in the same way that
`ilt` is called.

### References

See the source code for references.
