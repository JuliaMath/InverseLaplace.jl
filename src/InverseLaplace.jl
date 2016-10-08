module InverseLaplace

export Talbot, GWR

# Weeks
export  Weeks, WeeksErr, optimize, opteval

export ILtPair, abserr, iltpair_power

export talbot, gwr

abstract AbstractILt

type ILt{T<:Base.Callable, V<:Base.Callable} <: AbstractILt
    func::T
    iltfunc::V
    Nterms::Int
end

include("fixed_talbot.jl")
include("gwr.jl")
include("weeks.jl")
include("test.jl")

# Make this the default
ilt(args...) = talbot(args...)

# I think the vectorize macros already do this.
# And we should use the comprehensions anyway. They are faster for some reason.
# OTOH, ITL evaluation is relatively very slow.
for f in (:talbot, :gwr)
    @eval function $(f)(func, ta::AbstractArray, args...)
        length(ta) == 0 && return similar(ta)
        a1 = $(f)(func, ta[1], args...)
        a = similar(ta,typeof(a1))
        for (i,t) in enumerate(ta)
            a[i] = $(f)(func, t, args...)
        end
        a
    end
end

# Objects of type t are callable and use ILt function f
for (t,f) in ((:Talbot,:talbot), (:GWR, :gwr))
   @eval (ailt::$t)(t) = ($f)(ailt.func, t, ailt.Nterms)
   @eval (ailt::$t)(t,N) = ($f)(ailt.func, t, N)
end

end # module
