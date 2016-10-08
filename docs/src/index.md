# InverseLaplace
### Numerical inverse Laplace transform

Linux, OSX: [![Build Status](https://travis-ci.org/jlapeyre/InverseLaplace.jl.svg)](https://travis-ci.org/jlapeyre/InverseLaplace.jl)
&nbsp;
Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/jlapeyre/InverseLaplace.jl?branch=master&svg=true)](https://ci.appveyor.com/project/jlapeyre/inverselaplace-jl)
&nbsp; &nbsp; &nbsp;
[![Coverage Status](https://coveralls.io/repos/github/jlapeyre/InverseLaplace.jl/badge.svg?branch=master)](https://coveralls.io/github/jlapeyre/InverseLaplace.jl?branch=master)
[![codecov](https://codecov.io/gh/jlapeyre/InverseLaplace.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlapeyre/InverseLaplace.jl)

```@contents
```

## Inverse Laplace transform types

Constructing these types returns a callable object that evaluates the inverse transform at specified points.

```@docs
ILt
Weeks
WeeksErr
```

## Setting parameters

The inverse Laplace tranform routines are not black boxes. They are prone to instability and can give inaccurate or
wrong results. There are some parameters you can set to try to minimize these problems.

```@docs
setNterms
optimize
opteval
setparameters
```

## Investigating performance

```@docs
ILtPair
abserr
```

## Lower-level interface

Some of the lower-level routines can be called directly, without constructing types defined in `InverseLaplace`.

```@docs
ilt
talbot
gwr
InverseLaplace.talbotarr
```

## References

J.A.C Weideman, *Algorithms for Parameter Selection in the Weeks Method for Inverting the Laplace Transform,
SIAM Journal on Scientific Computing*, Vol. 21, pp. 111-128 **(1999)**


Abate, J. and Valkó, P.P., *Multi-precision Laplace transform inversion
International Journal for Numerical Methods in Engineering*, Vol. 60 (Iss. 5-7) pp 979–993 **(2004)**

Valkó, P.P. and Abate, J.,
*Comparison of Sequence Accelerators for the Gaver Method of Numerical Laplace Transform Inversion*,
Computers and Mathematics with Application,  Vol. 48 (Iss.3-40) pp. 629-636 **(2004)**
