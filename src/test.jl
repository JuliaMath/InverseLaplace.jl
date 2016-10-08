type TransformPair
    ft
    fs
end

type ILtPair{T} <: AbstractILt
    ilt::T
    ft
end

function iltpair_power(n)
    fs = s -> s^(- n - 1) * gamma(1+n)
    ft = t -> t^n
    TransformPair(ft,fs)
end

abserr(p::ILtPair, t) = abs(p.ilt(t) - p.ft(t))

# Create ILtPair with ilt type
for ilttype in (:Talbot, :GWR, :Weeks)
    @eval $(ilttype)(p::TransformPair,args...) = ILtPair($(ilttype)(p.fs,args...), p.ft)
end

for f in (:optimize,)
    @eval $(f)(p::ILtPair,args...) = $(f)(p.ilt, args...)
end
