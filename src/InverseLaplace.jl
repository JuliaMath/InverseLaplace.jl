module InverseLaplace

export Talbot, GWR, Weeks

#export ilt, gwr

abstract AbstractILT


include("fixed_talbot.jl")
include("gwr.jl")
include("weeks.jl")

# I think the vectorize macros already do this.
# And we should use the comprehensions anyway. They are faster for some reason.
for f in (:ilt, :gwr)
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

for (t,f) in ((:Talbot,:ilt), (:GWR, :gwr))
   @eval (ailt::$t)(t) = ($f)(ailt.func, t, ailt.Nterms)
   @eval (ailt::$t)(t,N) = ($f)(ailt.func, t, N)
end


end # module
