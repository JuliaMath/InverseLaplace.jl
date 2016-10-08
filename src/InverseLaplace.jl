module InverseLaplace

export ILt, setNterms

export  Weeks, WeeksErr, optimize, opteval, setparameters

export ILtPair, abserr, iltpair_power

export ilt, talbot, gwr

abstract AbstractILt

type ILt{T<:Base.Callable, V<:Base.Callable} <: AbstractILt
    func::T
    iltfunc::V
    Nterms::Int
end

"""
    itrans = ILt(func, iltfunc, Nterms=32)

return an object that estimates the inverse Laplace transform of
the function `func` using the algorithm implemented by function `iltfunc`.
`itrans(t)` estimates the inverse transform for argument `t`.  The
accuracy of the estimates depends strongly on the choice of `iltfunc`, `t`, `Nterms`,
and the precision of the data type of the argument to `func`. The default value
of `32` may give extremely inaccurate estimates.

`iltfunc` may be either `talbot` or `gwr`.
"""
ILt(func,iltfunc) = ILt(func, iltfunc, 32)

# Make this the default
"""
    ilt(func::Function, t::AbstractFloat, M::Integer=32)

`ilt` is an alias for the default inverse Laplace transform method `talbot`.
"""
ilt(args...) = talbot(args...)
ILt(func) = ILt(func, talbot, 32)


"""
    setNterms{T<:AbstractILt}(ailt::T, Nterms::Integer)

set the number of terms used in the inverse Laplace tranform `itrans`. If
`ailt` stores internal data, it will be recomputed, so that subsequent
calls `ailt(t)` reflect the new value of `Nterms`.
"""
setNterms(ailt::ILt, N::Integer) = (ailt.Nterms = N)

(ailt::ILt)(t) = ailt.iltfunc(ailt.func, t, ailt.Nterms)
(ailt::ILt)(t,N) = ailt.iltfunc(ailt.func, t, N)

include("fixed_talbot.jl")
include("gwr.jl")
include("weeks.jl")
include("test.jl")

# Broadcasting currently does not work because the first arg is a function
# talbot.( s -> 1/s^2 , [1.0,2.0])
# ERROR: MethodError: no method matching size(::##9#10)
#
for f in (:talbot, :gwr)
    @eval $(f)(func, ta::AbstractArray, args...) =  [ $(f)(func, t, args...) for t in ta ]
end

end # module
