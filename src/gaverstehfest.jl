# Gaver-Stehfest method

# 1.
# Stehfest, Harald.
# Algorithm 368: Numerical inversion of Laplace transforms [D5].
# Communications of the ACM 13.1 (1970): 47-49.

# 2.
# Kuhlman, Kristopher L.
# Review of inverse Laplace transform algorithms for Laplace-space numerical approaches.
# Numerical Algorithms 63.2 (2013): 339-355.

"""
    stehfest_coeffs(N::Integer)

Computes the Stehfest coefficients (summation weights) Vk for N terms. Coefficients are computed using BigFloats.
Only depends on the approximation order (N) and the precision. N must be even but can be precomputed and used only once.
Coefficients get very large in magnitude and oscillate rapidly. Must be careful to avoid catastrophic cancellation for lower precision.
"""
function stehfest_coeffs(N::Integer)
    if isodd(N)
        N += 1
        @warn "N must be even... increment N += 1"
    end
    v = zeros(BigFloat, N)
    aux = zero(eltype(v))
    for k in 1:N
        for j in fld(k + 1, 2):minimum(Int, [k, div(N, 2)])
            aux = big(j)^(div(N, 2)) * factorial(big(2 * j))
            aux /= factorial(big(div(N, 2) - j)) * factorial(big(j)) * factorial(big(j - 1))
            aux /= factorial(big(k - j)) * factorial(big(2 * j - k))
            v[k] += aux
        end
        v[k] *= (-1)^(k + div(N, 2))
    end
    return v
end

"""
    gaverstehfest(func::Function, t::AbstractFloat, v = stehfest_coeffs(N))

Evaluate the inverse Laplace transform of `func` at the point `t` using the Gaver-Stehfest algorithm.
`v` are the stehfest coefficients computed with `stehfest_coeffs` that only depend on number of terms. 
N (which must be even) defaults to 36 should depend on the precision used and desired accuracy.
The precision of the stehfest coefficients is used in the computation. 

In double precision (Float64), N should be <= 18 to provide best accuracy.
Usually only moderate accuracy can be achieved in double precision ~rtol ≈ 1e-4.
For higher accuracy, BigFloats should be used with N = 36.

Increasing precision should be accompanied by an increase in the number of coefficients used.
Increasing precision without increasing number of coefficients will not yield better accuracy. The inverse is generally true as well.

This method is not robust to oscillating F(t) and must be smooth.

# Examples
```
julia> F(s) = 1 / (s + 1) # where analytical inverse is f(t) = exp(-t)

julia> InverseLaplace.gaverstehfest(F, 2.0) # computes with default 36 coefficients
0.1353352832366128315426471959891170863784520272858265151734040491181311503059561

julia> InverseLaplace.gaverstehfest(F, 2.0, v = stehfest_coeffs(20)) # computes with custom number
0.1353353114073885136885007645878189977740624364169818512599342310325432753817199

# to calculate in double precision convert stehfest coefficients to Float64
julia> InverseLaplace.gaverstehfest(F, 2.0, v = convert(Vector{Float64}, stehfest_coeffs(18)))
0.13533650985980258
```
"""
function gaverstehfest(func::Function, t::AbstractFloat; v = stehfest_coeffs(36))
    N = length(v)
    a = zero(eltype(v))
    for k in 1:N
        bk = convert(eltype(v), k)
        a += v[k] * func(bk * log(convert(eltype(a), 2)) / t)
    end
    return a * log(convert(eltype(a), 2)) / t
end

"""
    gaverstehfest(func::Function, t::AbstractArray, v = stehfest_coeffs(N))

Evaluate the inverse Laplace transform of `func` over an array of `t` using the Gaver-Stehfest algorithm.
Computes coefficients once and calculates f(t) across available threads.

# Example
```
julia> gaverstehfest(s -> 1/(s + 1), 2.0:4.0)
0.1353355048178631463198819259249043857320738467297190227379428405365703748720082
0.04978728177016550841951683309410878070827215175940039571536203126773133604403101
0.01831383641619355549133471411401755204406266755309342902744499924574072011058623
```
"""
function gaverstehfest(f::Function, t::AbstractArray; v = stehfest_coeffs(18))
    N = length(v)
    a = zeros(eltype(v), length(t))
    Threads.@threads for ind in eachindex(t)
        a[ind] = gaverstehfest(f, t[ind], v = v)
    end
    return a
end
