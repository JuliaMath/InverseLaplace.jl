module InverseLaplace

export ILt, setNterms

#export Talbot, GWR, ILt, setNterms

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

ILt(func,iltfunc) = ILt(func, iltfunc, 32)

# Make this the default
ilt(args...) = talbot(args...)
ILt(func) = ILt(func, talbot, 32)

setNterms(ailt::ILt, N::Integer) = (ailt.Nterms = N)

(ailt::ILt)(t) = ailt.iltfunc(ailt.func, t, ailt.Nterms)
(ailt::ILt)(t,N) = ailt.iltfunc(ailt.func, t, N)

include("fixed_talbot.jl")
include("gwr.jl")
include("weeks.jl")
include("test.jl")

# I think the vectorize macros already do this.
# And we should use the comprehensions anyway. They are faster for some reason.
# OTOH, ITL evaluation is relatively very slow.
#
# Broadcasting currently does not work because the first arg is a function
# talbot.( s -> 1/s^2 , [1.0,2.0])
# ERROR: MethodError: no method matching size(::##9#10)
#
for f in (:talbot, :gwr)
    @eval $(f)(func, ta::AbstractArray, args...) =  [ $(f)(func, t, args...) for t in ta ]
end

#        length(ta) == 0 && return similar(ta)
        # a1 = $(f)(func, ta[1], args...)
        # a = similar(ta,typeof(a1))
        # for (i,t) in enumerate(ta)
        #     a[i] = $(f)(func, t, args...)
        # end
        # a
#    end

# phased out
# Objects of type t are callable and use ILt function f
# for (t,f) in ((:Talbot,:talbot), (:GWR, :gwr))
#    @eval (ailt::$t)(t) = ($f)(ailt.func, t, ailt.Nterms)
#    @eval (ailt::$t)(t,N) = ($f)(ailt.func, t, N)
# end

end # module
