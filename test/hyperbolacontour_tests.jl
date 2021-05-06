module hyperbolaTests

using Test
using InverseLaplace: hyperbola, hyper_fixed

#### F(s) = 1/sqrt(s), f(t) = 1/sqrt(π*t)

t = 1.0
@test isapprox(hyperbola(s -> 1/sqrt(s), t), 1/sqrt(pi*t))

t = 5.0
@test isapprox(hyperbola(s -> 1/sqrt(s), t), 1/sqrt(pi*t))

t = 1.0:10.0
@test isapprox(hyperbola(s -> 1/sqrt(s), t), 1 ./ sqrt.(pi.*t))
@test isapprox(hyper_fixed(s -> 1/sqrt(s), t), 1 ./ sqrt.(pi.*t))



#### F(s) = 1/s^4, f(t) = t^3 / 6

# test for single time point
t = 1.0
@test isapprox(hyperbola(s -> 1/s^4, t), t^3 / 6)

t = 5.0
@test isapprox(hyperbola(s -> 1/s^4, t), t^3 / 6)

# test across t array
t = 1.0:10.0
@test isapprox(hyperbola(s -> 1/s^4, t), t.^3 ./ 6)
@test isapprox(hyper_fixed(s -> 1/s^4, t, N = 36), t.^3 ./ 6) # note the increase of N to 36 so test can pass


#### F(s) = 1/(s+1), f(t) = exp(-t)

# test for single time point
t = 1.0
@test isapprox(hyperbola(s -> 1/(s+1), t), exp(-t))

t = 5.0
@test isapprox(hyperbola(s -> 1/(s+1), t), exp(-t))

# test across t array
t = 1.0:10.0
@test isapprox(hyperbola(s -> 1/(s+1), t), exp.(-t))
@test isapprox(hyper_fixed(s -> 1/(s+1), t), exp.(-t)) 

end
