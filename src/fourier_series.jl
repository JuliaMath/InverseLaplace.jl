# Abate, J., Choudhury, G.L., Whitt, W.
# (2000) An Introduction to Numerical Transform Inversion and Its Application to Probability Models. 
# Computational Probability. International Series in Operations Research & Management Science, vol 24
# The Fourier-Series Method for Laplace Transforms

# Implementation by Ridout, M.S.
# (2009) Generating random numbers from a distribution specified by its Laplace transform. Statistics and Computing, 19, 439-450

function _fourier_series_preprocess(m, L, A, nburn)
    nterms = nburn + m * L
    seqbtL = nburn:L:nterms
    y = π * im * collect(1:nterms) / L
    expy = exp.(y)
    A2L = 0.5 * A / L
    expxt = exp(A2L) / L

    coef = Vector{Float64}(undef, m+1)
    @inbounds for k in 0:m
        coef[k+1] = binomial(m, k) / 2^m
    end

    return seqbtL, y, expy, A2L, expxt, coef
end

function _fourier_series_compute(func, t, seqbtL, y, expy, A2L, expxt, coef)
    x = A2L / t
    z = x .+ y ./ t
    ltx = func(x)
    ltzexpy = func.(z) .* expy
    par_sum = 0.5 * real(ltx) .+ cumsum(real(ltzexpy))

    expxt * sum(coef .* par_sum[seqbtL]) / t
end

"""
    fourier_series(func, t, m=11, L=1, A=19, nburn=38)

Evaluate the inverse Laplace transform of `func` at the point `t`.

# Example

```jldoctest
julia> InverseLaplace.fourier_series( s -> 1/s^3, 3.0)
4.499985907607361
```

!!! note

    This function evaluates `func` only for real arguments.

"""
function fourier_series(func, t, m=11, L=1, A=19, nburn=38)
    seqbtL, y, expy, A2L, expxt, coef = _fourier_series_preprocess(m, L, A, nburn)

    _fourier_series_compute(func, t, seqbtL, y, expy, A2L, expxt, coef)
end