
<a id='InverseLaplace-1'></a>

# InverseLaplace


<a id='Numerical-inverse-Laplace-transform-1'></a>

### Numerical inverse Laplace transform


Linux, OSX: [![Build Status](https://travis-ci.org/jlapeyre/InverseLaplace.jl.svg)](https://travis-ci.org/jlapeyre/InverseLaplace.jl) &nbsp; Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/jlapeyre/InverseLaplace.jl?branch=master&svg=true)](https://ci.appveyor.com/project/jlapeyre/inverselaplace-jl) &nbsp; &nbsp; &nbsp; [![Coverage Status](https://coveralls.io/repos/github/jlapeyre/InverseLaplace.jl/badge.svg?branch=master)](https://coveralls.io/github/jlapeyre/InverseLaplace.jl?branch=master) [![codecov](https://codecov.io/gh/jlapeyre/InverseLaplace.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jlapeyre/InverseLaplace.jl)

- [InverseLaplace](index.md#InverseLaplace-1)
    - [Inverse Laplace transform types](index.md#Inverse-Laplace-transform-types-1)
    - [Setting parameters](index.md#Setting-parameters-1)
    - [Investigating performance](index.md#Investigating-performance-1)
    - [Lower-level interface](index.md#Lower-level-interface-1)
    - [References](index.md#References-1)


<a id='Inverse-Laplace-transform-types-1'></a>

## Inverse Laplace transform types


Constructing these types returns a callable object that evaluates the inverse transform at specified points.

<a id='InverseLaplace.ILt' href='#InverseLaplace.ILt'>#</a>
**`InverseLaplace.ILt`** &mdash; *Type*.



```
itrans = ILt(func, iltfunc, Nterms=32)
```

return an object that estimates the inverse Laplace transform of the function `func` using the algorithm implemented by function `iltfunc`. `itrans(t)` estimates the inverse transform for argument `t`.  The accuracy of the estimates depends strongly on the choice of `iltfunc`, `t`, `Nterms`, and the precision of the data type of the argument to `func`. The default value of `32` may give extremely inaccurate estimates.

`iltfunc` may be either `talbot` or `gwr`.


<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/InverseLaplace.jl#L19-L30' class='documenter-source'>source</a><br>

<a id='InverseLaplace.Weeks' href='#InverseLaplace.Weeks'>#</a>
**`InverseLaplace.Weeks`** &mdash; *Type*.



w::Weeks = Weeks(func::Function, Nterms::Integer=64, sigma=1.0, b=1.0)

return `w`, which estimates the inverse Laplace transform of `func` with the Weeks algorithm. `w(t)` evaluates the transform at `t`. The accuracy depends on the choice of `sigma` and `b`, with the optimal choices depending on `t`.

The call to `Weeks` that creates `w` is expensive relative to evaluation via `w(t)`.

**Example**

Compute the inverse transform of the transform of `cos` at argument `pi/2`.

```
julia> ft = Weeks(s -> s/(s^2+1), 80);

julia> ft(pi/2)
0.0
```


<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/weeks.jl#L28-L46' class='documenter-source'>source</a><br>

<a id='InverseLaplace.WeeksErr' href='#InverseLaplace.WeeksErr'>#</a>
**`InverseLaplace.WeeksErr`** &mdash; *Type*.



w::WeeksErr = WeeksErr(func::Function, Nterms::Integer=64, sigma=1.0, b=1.0)

return `w`, which estimates the inverse Laplace transform of `func` via the Weeks algorithm. `w(t)` returns a tuple containing the inverse transform at `t` and an error estimate. The accuracy of the inversion depends on the choice of `sigma` and `b`.

**Example**

Compute the inverse transform of the transform of `cos` at argument `pi/2` using `80` terms and an error estimate.

```
julia> ft = Weeks(s -> s/(s^2+1), 80);

julia> ft(pi/2)
(0.0,3.0872097665938698e-15)
```

This estimate is more accurate than `cos(pi/2)`.

```
julia> ft(pi/2)[1] - cos(pi/2)
-6.123233995736766e-17

julia> ft(pi/2)[1] - 0.0         # exact value
0.0

julia> ft(pi/2)[1] - cospi(1/2)  # cospi is more accurate
0.0
```


<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/weeks.jl#L162-L190' class='documenter-source'>source</a><br>


<a id='Setting-parameters-1'></a>

## Setting parameters


The inverse Laplace tranform routines are not black boxes. They are prone to instability and can give inaccurate or wrong results. There are some parameters you can set to try to minimize these problems.

<a id='InverseLaplace.setNterms' href='#InverseLaplace.setNterms'>#</a>
**`InverseLaplace.setNterms`** &mdash; *Function*.



```
setNterms{T<:AbstractILt}(ailt::T, Nterms::Integer)
```

set the number of terms used in the inverse Laplace tranform `itrans`. If `ailt` stores internal data, it will be recomputed, so that subsequent calls `ailt(t)` reflect the new value of `Nterms`.


<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/InverseLaplace.jl#L43-L49' class='documenter-source'>source</a><br>

<a id='InverseLaplace.optimize' href='#InverseLaplace.optimize'>#</a>
**`InverseLaplace.optimize`** &mdash; *Function*.



```
optimize{T<:AbstractWeeks}(w::T, t, Nterms)
```

optimize the parameters of the inverse Laplace transform `w` at the argument `t`. If `Nterms` is ommitted, the current value of `w.Nterms` is retained.

The accuracy of the Weeks algorithm depends strongly on `t`. For some ranges of `t`, the accuracy is relatively insensitive to the parameters. For other values of `t`, even using optimized parameters results in estimates of the inverse transform that are extremely inaccurate.

`optimize` is expensive in CPU time and allocation, it performs nested single-parameter optimization over two parameterss.


<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/weeks.jl#L65-L79' class='documenter-source'>source</a><br>

<a id='InverseLaplace.opteval' href='#InverseLaplace.opteval'>#</a>
**`InverseLaplace.opteval`** &mdash; *Function*.



```
opteval{T<:AbstractWeeks}(w::T, t, Nterms)
```

estimate an inverse Laplace transform at argument `t` using `w` after optimizing the parameters for `t`. If `Nterms` is omitted, then the current value of `w.Nterms` is used.

Use `Weeks` or `WeeksErr` to create `w`.


<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/weeks.jl#L90-L98' class='documenter-source'>source</a><br>

<a id='InverseLaplace.setparameters' href='#InverseLaplace.setparameters'>#</a>
**`InverseLaplace.setparameters`** &mdash; *Function*.



```
setparameters{T<:AbstractWeeks}(w::T, sigma, b, Nterms)
```

Set the parameters for the inverse Laplace transform object `w` and recompute the internal data. Subsequent calls `w(t)` will use these parameters. If `Nterms` or both `Nterms` and `b` are omitted, then their current values are retained.


<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/weeks.jl#L106-L112' class='documenter-source'>source</a><br>


<a id='Investigating-performance-1'></a>

## Investigating performance

<a id='InverseLaplace.ILtPair' href='#InverseLaplace.ILtPair'>#</a>
**`InverseLaplace.ILtPair`** &mdash; *Type*.



```
p = ILtPair(ilt::AbstractILt, ft::Function)
```

return an object of type `ILtPair` that associates `ilt` the inverse Laplace transform of a function with it "exact" numerical inverse `ft`. Calling `abserr(p,t)` returns the absolute error between the inverse transform and the exact value.        

**Example**

This example compares the inversion using the Weeks algorithm of the Laplace transform of `cos(t)` to its exact value at `t=1.0`.

```
julia> p = ILtPair( Weeks( s -> s/(1+s^2) ), cos);
julia> abserr(p,1.0)

0.0
```


<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/test.jl#L6-L24' class='documenter-source'>source</a><br>

<a id='InverseLaplace.abserr' href='#InverseLaplace.abserr'>#</a>
**`InverseLaplace.abserr`** &mdash; *Function*.



```
abserr(p::ILtPair, t)
```

Compute the absolute error between the estimated inverse Laplace transform and "exact" numerical solution contained in `p` at the point `t`. See `ILtPair`.


<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/test.jl#L36-L41' class='documenter-source'>source</a><br>


<a id='Lower-level-interface-1'></a>

## Lower-level interface


Some of the lower-level routines can be called directly, without constructing types defined in `InverseLaplace`.

<a id='InverseLaplace.ilt' href='#InverseLaplace.ilt'>#</a>
**`InverseLaplace.ilt`** &mdash; *Function*.



```
ilt(func::Function, t::AbstractFloat, M::Integer=32)
```

`ilt` is an alias for the default inverse Laplace transform method `talbot`.


<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/InverseLaplace.jl#L34-L38' class='documenter-source'>source</a><br>

<a id='InverseLaplace.talbot' href='#InverseLaplace.talbot'>#</a>
**`InverseLaplace.talbot`** &mdash; *Function*.



```
talbot(func::Function, t::AbstractFloat, M::Integer=32)
```

Evaluate the inverse Laplace transform of `func` at the point `t`. Use `M` terms in the algorithm. For `typeof(t)` is `Float64`, the default for `M` is `32`. For `BigFloat` the default is `64`.

If `BigFloat` precision is larger than default, try increasing `M`. `talbot is vectorized over`t`.

**Example**

```jlcon
julia> talbot( s -> 1/s^3,  3)
4.50000000000153
```

!!! note
    This function uses the fixed Talbot method. It evaluates `func` for complex arguments.



<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/fixed_talbot.jl#L13-L30' class='documenter-source'>source</a><br>

<a id='InverseLaplace.gwr' href='#InverseLaplace.gwr'>#</a>
**`InverseLaplace.gwr`** &mdash; *Function*.



```
gwr(func::Function, t::AbstractFloat, M::Integer=16)
```

Evaluate the inverse Laplace transform of `func` at the point `t`. Use `M` terms in the algorithm. For `typeof(t)` is `Float64`, the default for `M` is `16`. For `BigFloat` the default is `64`.

If `BigFloat` precision is larger than default, try increasing `M`.

**Example**

```jlcon
julia> gwr( s -> 1/s^3,  3.0)
4.499985907607361
```

!!! note
    This function uses the Gaver-Wynn rho method. It evaluates `func` only for real arguments.



<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/gwr.jl#L13-L33' class='documenter-source'>source</a><br>

<a id='InverseLaplace.talbotarr' href='#InverseLaplace.talbotarr'>#</a>
**`InverseLaplace.talbotarr`** &mdash; *Function*.



```
talbotarr{T}(func, ta::AbstractArray{T}, M)
```

inverse Laplace transform vectorized over `ta`. Each evaluation of `func(s)` is used for all elements of `ta`. This may be faster than a vectorized application of `talbot`, but is in general, less accurate. `talbotarr` uses the "fixed" Talbot method.


<a target='_blank' href='https://github.com/jlapeyre/InverseLaplace.jl/tree/d31f16b19c3924b0a4f96ed19f9fa5c63df52c8d/src/fixed_talbot.jl#L56-L63' class='documenter-source'>source</a><br>


<a id='References-1'></a>

## References


J.A.C Weideman, *Algorithms for Parameter Selection in the Weeks Method for Inverting the Laplace Transform, SIAM Journal on Scientific Computing*, Vol. 21, pp. 111-128 **(1999)**


Abate, J. and Valkó, P.P., *Multi-precision Laplace transform inversion International Journal for Numerical Methods in Engineering*, Vol. 60 (Iss. 5-7) pp 979–993 **(2004)**


Valkó, P.P. and Abate, J., *Comparison of Sequence Accelerators for the Gaver Method of Numerical Laplace Transform Inversion*, Computers and Mathematics with Application,  Vol. 48 (Iss.3-40) pp. 629-636 **(2004)**

