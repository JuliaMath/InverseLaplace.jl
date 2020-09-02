# Fixed Talbot method
#
# Abate, J. and Valkó, P.P.
# Multi-precision Laplace transform inversion
# International Journal for Numerical Methods in Engineering, Vol. 60 (Iss. 5-7)  2004  pp 979–993

"""
    talbot(func::Function, t::AbstractFloat, M::Integer=talbot_default_num_terms)

Evaluate the inverse Laplace transform of `func` at the point `t`. Use `M` terms in the algorithm.
For `typeof(t)` is `Float64`, the default for `M` is `32`. For `BigFloat` the default is `64`.

If `BigFloat` precision is larger than default, try increasing `M`.

# Example

```jldoctest
julia> InverseLaplace.talbot(s -> 1 / s^3, 3)
4.50000000000153
```

!!! note
    This function uses the fixed Talbot method. It evaluates `func` for complex arguments.
"""
function talbot(func, t, M)
    bM = convert(typeof(t), M)
    r = (2 * bM) / (5 * t)
    term = (1//2) * exp(r * t) * func(r)
    for i in 1:M-1
        theta = i * (pi/bM)
        s = r * theta * (complex(cot(theta), one(theta)))
        sigma = theta + (theta * cot(theta) - 1) * cot(theta)
        term += real(exp(t * s) * complex(one(t), sigma) * func(s))
    end
    return term * 2 / (5 * t)
end

talbot(func, t) = talbot(func, t, talbot_default_num_terms)
const talbot_BigFloat_default_num_terms = 64
talbot(func, t::BigFloat) = talbot(func, t, talbot_BigFloat_default_num_terms)
talbot(func, t::Integer, args...) = talbot(func, BigFloat(t), args...)
# Hmm, at some point, one of these routines actually gave a Rational result. Don't recall how.
# But, it can't be talbot.
talbot(func, t::Rational, args...) = talbot(func, BigFloat(t), args...)

# Operate on an array of values of t. A single function evaluation
# f(s) is used for all t
# This gives more inaccurate results the further values of
# t are from tmax
"""
    talbotarr(func, ta::AbstractArray, M)

Compute the inverse Laplace transform for each element in `ta`. Each evaluation
of `func(s)` is used for all elements of `ta`. This may be faster than a
broadcast application of `talbot` (i.e. `talbot.(...`) , but is in general, less
accurate.  `talbotarr` uses the "fixed" Talbot method.
"""
function talbotarr(func, t::AbstractArray, M)
    tt = typeof(t[1])
    bM = convert(tt, M)
    terms = similar(t)
    tmax = maximum(t)
    r = (2 * bM) / (5 * tmax)
    fr = (1//2) * func(r)
    terms .= exp.(r .* t) .* fr
    for i in 1:M-1
        theta = i * (pi/bM)
        s = r * theta * (complex(cot(theta), one(theta)))
        sigma = theta + (theta * cot(theta) - 1) * cot(theta)
        fs = complex(one(tt), sigma) * func(s)
        for j in 1:length(terms)
            terms[j] += real(exp(t[j] * s) * fs)
        end
    end
    terms .= terms .* 2 / (5 * tmax)
    return terms
end

talbotarr(func, t::AbstractArray) = talbotarr(func, t, talbot_default_num_terms)
talbotarr(func, t::Vector{BigFloat}) = talbotarr(func, t, talbot_BigFloat_default_num_terms)

#  LocalWords:  talbotarr func talbot
