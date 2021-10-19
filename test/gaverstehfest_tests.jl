module postwidderTests

using Test
using InverseLaplace: gaverstehfest, stehfest_coeffs

# Compare against results reported in Algorithm 368: Numerical inversion of Laplace transforms [D5].

#### F(s) = 1/sqrt(s), f(t) = 1/sqrt(π*t)

# test for single time point

## test in double precision Float64
t = 1.0
@test isapprox(gaverstehfest(s -> 1/sqrt(s), t, v = convert(Vector{Float64}, stehfest_coeffs(18))), 1/sqrt(pi*t), rtol = 1e-5)
t = 5.0
@test isapprox(gaverstehfest(s -> 1/sqrt(s), t, v = convert(Vector{Float64}, stehfest_coeffs(18))), 1/sqrt(pi*t), rtol = 1e-5)

# test in the default BigFloat
t = 1.0
@test isapprox(gaverstehfest(s -> 1/sqrt(s), t, v = stehfest_coeffs(36)), 1/sqrt(pi*t), rtol = 1e-16)
t = 5.0
@test isapprox(gaverstehfest(s -> 1/sqrt(s), t, v = stehfest_coeffs(36)), 1/sqrt(pi*t), rtol = 1e-16)

# test across t array in BigFloat
t = 1.0:10.0
@test isapprox(gaverstehfest(s -> 1/sqrt(s), t, v = stehfest_coeffs(36)), 1 ./ sqrt.(pi.*t), rtol = 1e-12)


#### F(s) = 1/s^4, f(t) = t^3 / 6

# test for single time point
## tests in double precision
t = 1.0
@test isapprox(gaverstehfest(s -> 1/s^4, t, v = convert(Vector{Float64}, stehfest_coeffs(18))), t^3/6, rtol = 1e-4)
t = 5.0
@test isapprox(gaverstehfest(s -> 1/s^4, t, v = convert(Vector{Float64}, stehfest_coeffs(18))), t^3/6, rtol = 1e-4)

# test in BigFloat
t = 1.0
@test isapprox(gaverstehfest(s -> 1/s^4, t, v = stehfest_coeffs(36)), t^3/6, rtol = 1e-14)
t = 5.0
@test isapprox(gaverstehfest(s -> 1/s^4, t, v = stehfest_coeffs(36)), t^3/6, rtol = 1e-14)

# test across t array
t = 1.0:10.0
@test isapprox(gaverstehfest(s -> 1/s^4, t, v = stehfest_coeffs(36)), t.^3 ./6, rtol = 1e-12)

#### F(s) = 1/(s+1), f(t) = exp(-t)
## for larger t you need a higher amount of coeffs depending on accuracy

# test for single time point
## double precision
t = 1.0
@test isapprox(gaverstehfest(s -> 1/(s+1), t, v = convert(Vector{Float64}, stehfest_coeffs(18))), exp(-t), rtol = 1e-4)
t = 5.0
@test isapprox(gaverstehfest(s -> 1/(s+1), t, v = convert(Vector{Float64}, stehfest_coeffs(18))), exp(-t), rtol = 1e-3)

# BigFloat
t = 1.0
@test isapprox(gaverstehfest(s -> 1/(s+1), t, v = stehfest_coeffs(36)), exp(-t), rtol = 1e-12)
t = 5.0
@test isapprox(gaverstehfest(s -> 1/(s+1), t, v = stehfest_coeffs(36)), exp(-t), rtol = 1e-8)

# test across t array
t = 1.0:10.0
@test isapprox(gaverstehfest(s -> 1/(s+1), t, v = stehfest_coeffs(36)), exp.(-t), rtol = 1e-8)

end # module
