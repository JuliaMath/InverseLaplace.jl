# Gaver-Stehfest method "Post-Widder"

# 1.
# Stehfest, Harald.
# Algorithm 368: Numerical inversion of Laplace transforms [D5].
# Communications of the ACM 13.1 (1970): 47-49.

# 2.
# Kuhlman, Kristopher L. 
# Review of inverse Laplace transform algorithms for Laplace-space numerical approaches.
# Numerical Algorithms 63.2 (2013): 339-355.

# 3. 
# Al-Shuaibi, Abdulaziz. 
# Inversion of the laplace transform via post—widder formula.
# Integral Transforms and Special Functions 11.3 (2001): 225-232.


#= Compute Post-Widder coefficients Vk for N terms. Also referred to as Salzer summation weights or Stehfest coefficients.
Only depend on the approximation order (N) and the precision. M must be even but can be precomputed and used only once. 
Coefficients get very large in magnitude and oscillate rapidly. Must be careful to avoid catastrophic cancellation for lower precision. 
For  N < 18 double precision is usually ok but less accurate, utilizng defualt BigFloat precision and N = 36 usually provides sufficient accuracy.
=#
function _PWcoeffs(N::Integer)
    if isodd(N)
        N += 1
        @warn "N must be even... increment N += 1"
    end
    v = zeros(BigFloat, N)
    aux = zero(eltype(v))
    for k in 1:N
        for j in floor(Int, div(k + 1, 2)):minimum(Int, [k, div(N, 2)])
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
    postwid(func::Function, t::AbstractFloat, v::Array{AbstractFloat,1}=_PWcoeffs(N))

Evaluate the inverse Laplace transform of `func` at the point `t` using the Gaver-Stehfest "Post-Widder" algorithm. 
Vk coefficients only depend on the number of expansion terms so can be computed once. N (which must be even) defaults to 18. 
In general the accuracy in double precision is poor, so arbitrary precision is used throughout. 

Increasing precision should be accompanied by an increase in the number of coefficients used.
Increasing precision without increasing number of coefficients will not yield better accuracy. The inverse is generally true as well.

This method is not robust to oscillating F(t) and must be smooth. 
# Example

```jldoctest
julia> setprecision(53); InverseLaplace.postwid(s -> 1/(s + 1), 2.0)
0.1353356835639731
```
"""
function postwid(func::Function, t::AbstractFloat; v = _PWcoeffs(18))
    N = length(v)
    a = zero(eltype(v))
    for k in 1:N
        bk = convert(eltype(v), k)
        a += v[k] * func(bk * log(convert(eltype(a), 2)) / t)
    end
    return a * log(convert(eltype(a), 2)) / t
end
"""
    postwid(func::Function, t::AbstractArray, v::Array{AbstractFloat,1}=_PWcoeffs(N))

Evaluate the inverse Laplace transform of `func` over an array of `t` using the Gaver-Stehfest "Post-Widder" algorithm. 
Computes coefficients once and calculates f(t) across available threads.
# Example

```jldoctest
julia> setprecision(53); postwid(s -> 1/(s + 1), 2.0:4.0)
 0.1353356835639731
 0.049786688998390505
 0.01831327193414956
```
"""
function postwid(f::Function, t::AbstractArray; v = _PWcoeffs(18))
    N = length(v)
    a = zeros(eltype(v), length(t))
    Threads.@threads for ind in eachindex(t)
        a[ind] = postwid(f, t[ind], v = v)
    end
    return a
end
