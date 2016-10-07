# Abate, J. and Valkó, P.P.
# Multi-precision Laplace transform inversion
# International Journal for Numerical Methods in Engineering, Vol. 60 (Iss. 5-7)  2004  pp 979–993
# Fixed Talbot method

type Talbot <: AbstractILT
    func::Function
    Nterms::Int
end

Talbot(func::Function) = Talbot(func, 32)

"""
    ilt(func::Function, t::AbstractFloat, M::Integer=32)

Evaluate the inverse Laplace transform of `func` at the point `t`. Use `M` terms in the algorithm.
For `typeof(t)` is `Float64`, the default for `M` is `32`. For `BigFloat` the default is `64`.

If `BigFloat` precision is larger than default, try increasing `M`. `ilt is vectorized over `t`.

# Example

```jldoctest
julia> ilt( s -> 1/s^3,  3)
4.50000000000153
```

!!! note
    This function uses the fixed Talbot method.
    It evaluates `func` for complex arguments.
"""
function ilt(func, t, M)
    bM = convert(typeof(t),M)
    r = (2 * bM) / (5*t)
    term = (1//2) * exp(r*t) * func(r)
    for i in 1:M-1
        theta = i * (pi/bM)
        s = r*theta*(complex(cot(theta),one(theta)))
        sigma = theta + (theta*cot(theta)-1)*cot(theta)
        term += real(exp(t*s) * complex(one(t),sigma) * func(s))
    end
    return term * 2 / (5*t)
end

ilt(func,t) = ilt(func,t,32)
ilt(func,t::BigFloat) = ilt(func,t,64)
ilt(func,t::Integer,args...) = ilt(func,BigFloat(t),args...)
# Hmm, at some point, one of these routines actually gave a Rational result. Don't recall how.
# But, it can't be ilt.
ilt(func,t::Rational,args...) = ilt(func,BigFloat(t),args...)


# Operate on an array of values of t. A single function evaluation
# f(s) is used for all t
# This gives more inaccurate results the further values of
# t are from tmax

"""
    iltarr{T}(func, ta::AbstractArray{T}, M)

inverse laplace transform vectorized over `ta`. Each evaluation
of `func(s)` is used for all elements of `ta`. This may be faster
than a vectorized application of `ilt`, but is in general, less accurate.
`iltarr` uses the "fixed" Talbot method.
"""
function iltarr{T}(func, t::AbstractArray{T}, M)
    tt = typeof(t[1])
    bM = convert(tt,M)
    terms = similar(t)
    tmax = maximum(t)
    r = (2 * bM) / (5*tmax)
    fr = (1//2) * func(r)
    for j in 1:length(terms)
        terms[j] =  exp(r*t[j]) * fr
    end
    for i in 1:M-1
        theta = i * (pi/bM)
        s = r*theta*(complex(cot(theta),one(theta)))
        sigma = theta + (theta*cot(theta)-1)*cot(theta)
        fs = complex(one(tt),sigma) * func(s)
        for j in 1:length(terms)
            terms[j] += real(exp(t[j]*s) * fs)
        end
    end
    for j in 1:length(terms)
        terms[j] *= 2/(5*tmax)
    end
    return terms
end

iltarr(func,t::AbstractArray) = iltarr(func,t,32)
iltarr(func,t::AbstractArray{BigFloat}) = iltarr(func,t,64)
