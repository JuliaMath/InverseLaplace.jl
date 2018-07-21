# InverseLaplace

*Numerical inverse Laplace transform*

The source repository is [https://github.com/jlapeyre/InverseLaplace.jl](https://github.com/jlapeyre/InverseLaplace.jl).

This package provides three algorithms for numerically inverting Laplace transforms.
`InverseLaplace` v0.1.0 is the last version that supports Julia v0.6.
Optimization of the Weeks method is temporarily disabled for Julia v0.7.

## Contents

```@contents
```

## Index

```@index
```

## Inverse Laplace transform algorithms

Constructing these Julia types, corresponding to different numerical algorithms,
returns a callable object that evaluates the inverse Laplace transform at specified points.

```@docs
Talbot
ILT
GWR
Weeks
WeeksErr
```

## Setting parameters

The inverse Laplace tranform routines should not be treated as black boxes. They are
prone to instability and can give inaccurate or wrong results. There are some
parameters you can set to try to minimize these problems.

```@docs
InverseLaplace.setNterms
InverseLaplace.optimize
InverseLaplace.opteval
InverseLaplace.setparameters
```

## Analzying accuracy

```@docs
ILtPair
abserr
iltpair_power
TransformPair
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
