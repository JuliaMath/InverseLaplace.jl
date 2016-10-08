type TransformPair
    ft
    fs
end

"""
    p = ILtPair(ilt::AbstractILt, ft::Function)

return an object of type `ILtPair` that associates `ilt` the inverse Laplace transform of
a function with it "exact" numerical inverse `ft`. Calling `abserr(p,t)` returns the
absolute error between the inverse transform and the exact value.        

# Example

This example compares the inversion using the Weeks algorithm of the Laplace transform
of `cos(t)` to its exact value at `t=1.0`.
```
julia> p = ILtPair( Weeks( s -> s/(1+s^2) ), cos);
julia> abserr(p,1.0)

0.0
```

"""
type ILtPair{T} <: AbstractILt
    ilt::T
    ft
end

function iltpair_power(n)
    fs = s -> s^(- n - 1) * gamma(1+n)
    ft = t -> t^n
    TransformPair(ft,fs)
end

"""
    abserr(p::ILtPair, t)

Compute the absolute error between the estimated inverse Laplace transform and
"exact" numerical solution contained in `p` at the point `t`. See `ILtPair`.
"""
abserr(p::ILtPair, t) = abs(p.ilt(t) - p.ft(t))

# Create ILtPair with ilt type
for ilttype in (:Talbot, :GWR, :Weeks)
    @eval $(ilttype)(p::TransformPair,args...) = ILtPair($(ilttype)(p.fs,args...), p.ft)
end

for f in (:optimize,)
    @eval $(f)(p::ILtPair,args...) = $(f)(p.ilt, args...)
end
