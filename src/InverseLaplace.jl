__precompile__()

module InverseLaplace

import SpecialFunctions
# I have to export everything in order for Documenter.jl to find
# the strings. A few hours of work would probably be enough to solve the problem.
export ILt, setNterms
# export optimize, opteval # Broken at the moment
export Weeks, WeeksErr, setparameters
export Talbot, GWR, ILT, Fourier_series
export TransformPair, ILtPair, abserr, iltpair_power
export ilt, talbot, gwr, fourier_series
export gaverstehfest, stehfest_coeffs

abstract type AbstractILt end

"""
   ft::Talbot = Talbot(func::Function, Nterms::Integer=32)

Return `ft`, which estimates the inverse Laplace transform of `func` with
the fixed Talbot algorithm. `ft(t)` evaluates the transform at `t`.  You may
want to tune `Nterms` together with `setprecision(BigFloat, x)`.

# Example

Compute the inverse transform of the transform of `cos` at argument `pi/2`.
```julia-repl
julia> ft = Talbot(s -> s / (s^2 + 1), 80);

julia> ft(pi / 2)
-3.5366510684573195e-5
```

Note that given `Float64` input, the precision of the returned value may not be satisfactory.
```julia-repl
julia> Float64(ft(big(pi) / 2))
2.114425886215604e-49
```

!!! note
    This function uses the fixed Talbot method. It evaluates the Laplace transform
    function for complex arguments. The GWR method is, in general, less accurate and
    less stable, but does not evaluate the Laplace transform function for complex arguments.
"""
struct Talbot{T<:Base.Callable} <: AbstractILt
    Laplace_space_func::T  # Laplace space function
    Nterms::Int
end

const talbot_default_num_terms = 32
Talbot(Laplace_space_func::Base.Callable) = Talbot(Laplace_space_func, talbot_default_num_terms)

"""
   ft::GWR = GWR(func::Function, Nterms::Integer=16)

Return `ft`, which estimates the inverse Laplace transform of `func` with
the GWR algorithm. `ft(t)` evaluates the transform at `t`.

# Example

Compute the inverse transform of the transform of `cos` at argument `pi / 2`.
```
julia> ft = GWR(s -> s / (s^2 + 1), 16);

julia> ft(pi / 2)
-0.001
```
"""
struct GWR{T<:Base.Callable} <: AbstractILt
    Laplace_space_func::T  # Laplace space function
    Nterms::Int
end

const gwr_default_num_terms = 16
GWR(Laplace_space_func::Base.Callable) = GWR(Laplace_space_func, gwr_default_num_terms)

"""
   ft::Fourier_series = Fourier_series(func::Function, m=11, L=1, A=19, nburn=38)

Return `ft`, which estimates the inverse Laplace transform of `func` with
the fourier series algorithm. `ft(t)` evaluates the transform at `t`.

# Example

Compute the inverse transform of the transform of `cos` at argument `pi / 2`.
```
julia> ft = Fourier_series(s -> s / (s^2 + 1), 16);

julia> ft(pi / 2)
-0.001
```
"""
struct Fourier_series{T<:Base.Callable} <: AbstractILt
    Laplace_space_func::T  # Laplace space function
    m::Integer
    L::Integer
    A::Integer
    nburn::Integer
end

const fourier_series_default_m = 11
const fourier_series_default_L = 1
const fourier_series_default_A = 19
const fourier_series_default_nburn = 38
Fourier_series(Laplace_space_func::Base.Callable) = Fourier_series(
    Laplace_space_func,
    fourier_series_default_m,
    fourier_series_default_L,
    fourier_series_default_A,
    fourier_series_default_nburn
)


"""
    ILT(function, Nterms=32)

This is an alias for the default `Talbot()` method.
"""
ILT(args...) = Talbot(args...)

function Base.show(io::IO, ailt::AbstractILt)
    print(io, string(typeof(ailt)), "(Nterms=", ailt.Nterms, ')')
end

# TODO: get rid of this in favor of above.
# Rely on broadcasting, as well.
struct ILt{T<:Base.Callable, V<:Base.Callable} <: AbstractILt
    Laplace_space_func::T
    iltfunc::V
    Nterms::Int
end

"""
    itrans = ILt(func, iltfunc=talbot, Nterms=32)

 * deprecated. Use, ILT, Talbot, GWR, or Weeks instead *

Return an object that estimates the inverse Laplace transform of
the function `func` using the algorithm implemented by function `iltfunc`.
`itrans(t)` estimates the inverse transform for argument `t`.  The
accuracy of the estimates depends strongly on the choice of `iltfunc`, `t`, `Nterms`,
and the precision of the data type of the argument to `func`. The default value
of `32` may give extremely inaccurate estimates.

`iltfunc` may be either `talbot` or `gwr`.

## Example

```julia
julia> itr = ILt( s -> 1/(1 + s^2), talbot);

julia> itr([ pi/4, pi/2, 3*pi/4, -pi])
4-element Array{Float64,1}:
  0.707107
  1.0
  0.707107
 -3.66676e-12
```
"""
ILt(func, iltfunc) = ILt(func, iltfunc, talbot_default_num_terms)

# Make this the default
"""
    ilt(func::Function, t::AbstractFloat, M::Integer=32)

`ilt` is an alias for the default inverse Laplace transform method `talbot`.
"""
ilt(args...) = talbot(args...)
ILt(func) = ILt(func, talbot, talbot_default_num_terms)

"""
    setNterms(ailt::AbstractILt, Nterms::Integer)

set the number of terms used in the inverse Laplace tranform `ailt`. If
`ailt` stores internal data, it will be recomputed, so that subsequent
calls `ailt(t)` reflect the new value of `Nterms`.
"""
setNterms(ailt::AbstractILt, N::Integer) = (ailt.Nterms = N)

## Make the ILT data types callable.
(ailt::ILt)(t) = ailt.iltfunc(ailt.Laplace_space_func, t, ailt.Nterms)
(ailt::ILt)(t,N) = ailt.iltfunc(ailt.Laplace_space_func, t, N)

(ailt::Talbot)(t) = talbot(ailt.Laplace_space_func, t, ailt.Nterms)
(ailt::Talbot)(t, n::Integer) = talbot(ailt.Laplace_space_func, t, n)

(ailt::GWR)(t) = gwr(ailt.Laplace_space_func, t, ailt.Nterms)
(ailt::GWR)(t, n::Integer) = gwr(ailt.Laplace_space_func, t, n)

(ailt::Fourier_series)(t) = fourier_series(ailt.Laplace_space_func, t)
(ailt::Fourier_series)(t, m::Integer, L::Integer, A::Integer, nburn::Integer) = fourier_series(ailt.Laplace_space_func, t, m, L, A, nburn)

include("fixed_talbot.jl")
include("gwr.jl")
include("weeks.jl")
include("pairtest.jl")
include("gaverstehfest.jl")
include("fourier_series.jl")

end # module
