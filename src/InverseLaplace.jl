module InverseLaplace

export ilt, gwr

# Two methods to compute the inverse Laplace transform.

# Abate, J. and Valkó, P.P.
# Multi-precision Laplace transform inversion
# International Journal for Numerical Methods in Engineering, Vol. 60 (Iss. 5-7)  2004  pp 979–993
# Fixed Talbot method

"""
    ilt(func::Function, t::AbstractFloat, M::Integer=32)

Evaluate the inverse Laplace transform of `func` at the point `t`. Use `M` terms in the algorithm.
For `typeof(t)` is `Float64`, the default for `M` is `32`. For `BigFloat` the default is `64`.

If `BigFloat` precision is larger than default, try increasing `M`.

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
ilt(func,t::Integer) = ilt(func,BigFloat(t))
ilt(func,t::Rational) = ilt(func,BigFloat(t))


# Operate on an array of values of t. A single function evaluation
# f(s) is used for all t
# This gives more inaccurate results the further values of
# t are from tmax
function ilt{T}(func, t::AbstractArray{T}, M)
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

ilt(func,t::AbstractArray) = ilt(func,t,32)
ilt(func,t::AbstractArray{BigFloat}) = ilt(func,t,64)

# Valkó, P.P. and Abate, J.
# Comparison of Sequence Accelerators for the Gaver Method of Numerical Laplace Transform Inversion
# Computers and Mathematics with Application,  Vol. 48 (Iss.3-40) 2004 pp. 629-636
# Gaver Wynn rho method

"""
    gwr(func::Function, t::AbstractFloat, M::Integer=16)

Evaluate the inverse Laplace transform of `func` at the point `t`. Use `M` terms in the algorithm.
For `typeof(t)` is `Float64`, the default for `M` is `16`. For `BigFloat` the default is `64`.

If `BigFloat` precision is larger than default, try increasing `M`.

# Example

```jldoctest
julia> gwr( s -> 1/s^3,  3.0)
4.499985907607361
```

!!! note
    This function uses the Gaver-Wunn rho method.
    It evaluates `func` only for real arguments.
"""
function gwr(func, t, M)
    Dt = typeof(t)
    bM = convert(Dt,M)
    tau = log(convert(Dt,2))/t
    broken = false
    Fi = Array(Dt,2 * M)
    for i in 1: 2 * M
        Fi[i] = func(i * tau)
    end
    M1 = M
    G0 = zeros(Dt,M1+1)
    for n in 1:M
        sm = zero(Dt)
        bn = convert(Dt,n)
        for i in 0:n
            bi = convert(Dt,i)
            sm += binomial(big(n),big(i)) * (-1)^i * Fi[n+i]
        end
        G0[n] = tau * factorial(2*bn)/(factorial(bn)*factorial(bn-1)) * sm
    end
    Gm = zeros(Dt,M1+1)
    Gp = zeros(Dt,M1+1)
    best = G0[M1]
    for k in 0:M1-2
        for n in (M1-2-k):-1:0
            expr = G0[n+2] - G0[n+1]
            if expr == 0
                broken = true
                break
            end
            expr = Gm[n+2] + (k+1)/expr
            Gp[n+1] = expr
            if isodd(k) && n == M1 - 2 - k
                best = expr
            end
        end
        if broken break end
        for n in 0:(M1-k)
            Gm[n+1] = G0[n+1]
            G0[n+1] = Gp[n+1]
        end
    end
    best
end

gwr(func, t::Float64) = gwr(func,t,16)
gwr(func, t::BigFloat) = gwr(func,t,64)
gwr(func, t::Integer) = gwr(func,BigFloat(t))
gwrt(func, t::Real) = gwr(func,float(t))

end # module
