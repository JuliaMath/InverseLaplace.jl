### Implement Laplace transform along a hyperbola contour ###

# 1.
# Weideman, J., and L. Trefethen.
# Parabolic and hyperbolic contours for computing the Bromwich integral.
# Mathematics of Computation 76.259 (2007): 1341-1356.
# 2.
# Liemert, André, and Alwin Kienle.
# Application of the Laplace transform in time-domain optical spectroscopy and imaging.
# Journal of biomedical optics 20.11 (2015): 110502.

"""
    hyperbola(f::Function, t::AbstractFloat; N::Int = 16)

Evaluate the inverse Laplace transform of `f` at the point `t` by approximating the Bromwhich integral with a hyperbola contour.
A real valued time-domain signal is assumed with the integral being calculated by a computable series applying the midpoint rule.

The contour is defined by several parameters. μ controls the extreme nodes with a larger parameter moving the outlier nodes further into left half plane.
h is the uniform node spacing with N being the number of nodes on the contour. The parameters depend on both t and N such that the function must be computed N times for each t.
This is very inefficient especially if f is computational expensive, so this approach is recommended when only a few values of t are needed.
However, this method will provide more accurate results at the cost of computational time compared to fixed contour approaches.

The number of nodes, N, defaults to 16. The optimal N is case dependent and arbitrarily increasing N will not result in more accurate results. Smaller values of N will be faster.

Bromwhich contour approaches should only be applied for t > 0. We are aiming to optimize the convergence as N -> inf so our parameters are ~ N / t.
In other words we want the Bromwhich integral to converge rapidly, and this can be done by starting and ending our integration path in the left hand plane causing exp(st) to decay.
If t = 0, the exponential factor does not cause good convergence. If you want to evaluate for t = 0, you should try a small positive value like 1e-30.

See also: [`hyper_fixed`](@ref), [`talbot`](@ref), [`talbotarr`](@ref)

# Example

```jldoctest
julia> InverseLaplace.hyperbola(s -> 1/(s + 1), 2.0)
0.13533528323665164
julia> InverseLaplace.hyperbola.(s -> 1/(s + 1), 1.0:3.0)
3-element Vector{Float64}:
 0.36787944117147775
 0.13533528323665164
 0.04978706836793606
"""
function hyperbola(f::Function, t::T; N::Int = 16) where T
    a =  zero(Complex{T})
    h = T(1.081792140) / N
    for k in 0:N-1
        sk, dsk = s((k + T(0.5)) * h, N, t)
        a += f(sk) * exp(sk * t) * dsk
    end
    return imag(a) * h / π
end

# adaptive contour for each time point (can't precompute f(s) or sk)
# parameterized hyperbola contour
# parameters are given by ref [2] which are given in increased precision from calculations given in [1]
function s(θ, N, t::T) where T
    μ = T(4.492075287) * N / t
    ϕ = T(1.172104229)
    a = θ + im * ϕ
    s = μ + im * μ * sinh(a)
    ds = im * μ * cosh(a) # derivitive of hyperbola contour
    return s, ds
end

"""
    hyper_fixed(f::Function, t::AbstractVector; N::Int = 24)

Evaluate the inverse Laplace transform of `f` over a vector `t` by approximating the Bromwhich integral with a fixed hyperbola contour.
A real valued time-domain signal is assumed with the integral being calculated by a computable series applying the midpoint rule.

In contrast to `hyperbola` where the contour is dependent on 't', the contour is fixed over the entire vector of `t` values.
This function should be used when F(s) is computationally expensive to minimize the number of evaulations of F(s) or when the transform
needs to be applied over many values of `t`.

It operates under the assumption that f(t) is needed for many values of `t` over some interval `t ∈ (first(t):t[end])` and t>0.
The number of nodes `N` can be increased for accuracy, however there is an optimal N that minimizes accuracy where error increases afterwords.

This approach will in general be less accurate than `hyperbola` as we are using the same contour over an array of time values, instead of an optimized contour at each time.
Bromwhich contour approaches should only be applied for t > 0.

See also: [`hyper_fixed!`](@ref), [`hyperbola`](@ref), [`talbot`](@ref), [`talbotarr`](@ref)

# Example

```jldoctest
julia> InverseLaplace.hyper_fixed(s -> 1/(s + 1), 2.0:3.0)
2-element Vector{Float64}:
 0.13533528323660948
 0.0497870683678199
"""
hyper_fixed(f::Function, t::AbstractVector{T}; N::Int = 24) where T = hyper_fixed!(zeros(T, length(t)), f, t, N = N)

"""
    hyper_fixed!(out::AbstractVector{T}, f::Function, t::AbstractVector{T}; N::Int = 24)

Evalue the inverse Laplace transform of `f` over a vector of `t` using a fixed hyperbola contour.
The result is stored in `out` which must be distinct from `t` (they cannot alias each other) and be a vector of zeros.

See [`hyper_fixed`](@ref) for more details about the implementation.

# Example

```jldoctest
julia> t = 1.0:0.2:2.0;
julia> out = zeros(eltype(t), length(t));
julia> InverseLaplace.hyper_fixed!(out, s -> 1/(s + 1), t)
6-element Vector{Float64}:
 0.3678794411714416
 0.3011942119122033
 0.24659696394161384
 0.20189651799466898
 0.16529888822160033
 0.13533528323663216
"""
function hyper_fixed!(out::AbstractVector{T}, f::Function, t::AbstractVector{T}; N::Int = 24) where T
    length(out) == length(t) || throw(DimensionMismatch("out is not equal to length of t"))
    iszero(out) || throw(ArgumentError("out is not a vector of zeros"))

    μ, h = hyper_coef(N, t)

    for k in 0:N-1
        kh = (k + 1//2) * h
        sk = s_fixed(kh, μ)
        dsk = ds_fixed(kh, μ)
        a = f(sk) * dsk * h / π
        for ind in eachindex(t)
            out[ind] += imag(a * exp(sk * t[ind]))
        end
    end
    return out
end

# get the fixed integration components
function hyper_coef(N, t::AbstractVector{T}) where T
    ϕ = T(1.09)
    a = (T(π) - 2 * ϕ) * t[end] / first(t) + 4 * ϕ - T(π)
    a /= (4 * ϕ - π) * sin(ϕ)
    A = acosh(a)
    μ = (4 * T(π) * ϕ - T(π)^2) * N / t[end] / A
    h = A / N
    return μ, h
end

# compute the fixed hyperbola contour
s_fixed(θ, μ; ϕ = 1.09) = μ + im * μ * sinh(θ + im * ϕ)
ds_fixed(θ, μ; ϕ = 1.09) = im * μ * cosh(θ + im * ϕ)
