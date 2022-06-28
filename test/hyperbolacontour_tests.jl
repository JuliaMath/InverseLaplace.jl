module hyperbolaTests

using Test
using InverseLaplace: hyperbola, hyper_fixed

## test type stability

@test hyperbola(s -> 1/sqrt(s), 1.0f0) isa Float32
@test hyperbola(s -> 1/sqrt(s), 1.0) isa Float64
@test hyperbola(s -> 1/sqrt(s), big"1.0") isa BigFloat

@test hyper_fixed(s -> 1/sqrt(s), [1.0f0, 2.0f0])[1] isa Float32
@test hyper_fixed(s -> 1/sqrt(s), [1.0, 2.0])[1] isa Float64
@test hyper_fixed(s -> 1/sqrt(s), BigFloat.([1.0, 2.0]))[1] isa BigFloat

#### F(s) = 1/sqrt(s), f(t) = 1/sqrt(π*t)

t = 1.0
@test isapprox(hyperbola(s -> 1/sqrt(s), t), 1/sqrt(pi*t))

t = 5.0
@test isapprox(hyperbola(s -> 1/sqrt(s), t), 1/sqrt(pi*t))

t = 1.0:10.0
@test isapprox(hyperbola.(s -> 1/sqrt(s), t), 1 ./ sqrt.(pi.*t))
@test isapprox(hyper_fixed(s -> 1/sqrt(s), t), 1 ./ sqrt.(pi.*t))


#### F(s) = 1/s^4, f(t) = t^3 / 6

# test for single time point
t = 1.0
@test isapprox(hyperbola(s -> 1/s^4, t), t^3 / 6)

t = 5.0
@test isapprox(hyperbola(s -> 1/s^4, t), t^3 / 6)

# test across t array
t = 1.0:10.0
@test isapprox(hyperbola.(s -> 1/s^4, t), t.^3 ./ 6)
@test isapprox(hyper_fixed(s -> 1/s^4, t, N = 36), t.^3 ./ 6) #note the increase of N to 36 so test can pass


#### F(s) = 1/(s+1), f(t) = exp(-t)

# test for single time point
t = 1.0
@test isapprox(hyperbola(s -> 1/(s+1), t), exp(-t))

t = 5.0
@test isapprox(hyperbola(s -> 1/(s+1), t), exp(-t))

# test across t array
t = 1.0:10.0
@test isapprox(hyperbola.(s -> 1/(s+1), t), exp.(-t))
@test isapprox(hyper_fixed(s -> 1/(s+1), t), exp.(-t))

#### F(s) = 1/s - 1/(s + 1), f(t) = 1 - exp(-t)

# test for single time point
t = 1.0
@test isapprox(hyperbola(s -> 1/s - 1/(s + 1), t), 1 - exp(-t))

t = 5.0
@test isapprox(hyperbola(s -> 1/s - 1/(s + 1), t), 1 - exp(-t))

# test across t array
t = 1.0:10.0
@test isapprox(hyperbola.(s -> 1/s - 1/(s + 1), t), 1 .- exp.(-t))
@test isapprox(hyper_fixed(s -> 1/s - 1/(s + 1), t), 1 .- exp.(-t))

# evaluation at t = 0 requires inputting a positive value close to 0
t = 1e-30
@test isapprox(hyperbola(s -> 1/s - 1/(s + 1), t), 1 - exp(-0))

t = 1e-2:10.0
@test isapprox(hyperbola.(s -> 1/s - 1/(s + 1), t), 1 .- exp.(-t))
@test isapprox(hyper_fixed(s -> 1/s - 1/(s + 1), t, N = 56), 1 .- exp.(-t))

end
