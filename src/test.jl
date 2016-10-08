type TransformPair
    ft
    fs
end

type ILTPair{T} <: AbstractILt
    ilt::T
    ft
end

function iltpair_power(n)
    fs = s -> s^(- n - 1) * gamma(1+n)
    ft = t -> t^n
    TransformPair(ft,fs)
end

abserr(p::ILTPair, t) = abs(p.ilt(t) - p.ft(t))

# Create ILTPair with ilt type
for ilttype in (:Talbot, :GWR, :Weeks)
    @eval $(ilttype)(p::TransformPair,args...) = ILTPair($(ilttype)(p.fs,args...), p.ft)
end

