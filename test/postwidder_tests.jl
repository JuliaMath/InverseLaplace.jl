module postwidderTests

using Test
using InverseLaplace: postwid, _PWcoeffs


# Compare against results reported in Algorithm 368: Numerical inversion of Laplace transforms [D5].

#### F(s) = 1/sqrt(s), f(t) = 1/sqrt(π*t)

# test for single time point
setprecision(53)
t = 1.0
@test isapprox(postwid(s -> 1/sqrt(s), t, v = _PWcoeffs(18)), 1/sqrt(pi*t), rtol = 1e-4)
t = 5.0
@test isapprox(postwid(s -> 1/sqrt(s), t, v = _PWcoeffs(18)), 1/sqrt(pi*t), rtol = 1e-4)

setprecision(128)
t = 1.0
@test isapprox(postwid(s -> 1/sqrt(s), t, v = _PWcoeffs(36)), 1/sqrt(pi*t), rtol = 1e-12)
t = 5.0
@test isapprox(postwid(s -> 1/sqrt(s), t, v = _PWcoeffs(36)), 1/sqrt(pi*t), rtol = 1e-12)

# test across t array
t = 1.0:10.0

setprecision(53)
@test isapprox(postwid(s -> 1/sqrt(s), t, v = _PWcoeffs(18)), 1 ./ sqrt.(pi.*t), rtol = 1e-4)

setprecision(128)
@test isapprox(postwid(s -> 1/sqrt(s), t, v = _PWcoeffs(36)), 1 ./ sqrt.(pi.*t), rtol = 1e-12)


#### F(s) = 1/s^4, f(t) = t^3 / 6

# test for single time point
setprecision(53)
t = 1.0
@test isapprox(postwid(s -> 1/s^4, t, v = _PWcoeffs(18)), t^3/6, rtol = 1e-4)
t = 5.0
@test isapprox(postwid(s -> 1/s^4, t, v = _PWcoeffs(18)), t^3/6, rtol = 1e-4)

setprecision(128)
t = 1.0
@test isapprox(postwid(s -> 1/s^4, t, v = _PWcoeffs(36)), t^3/6, rtol = 1e-12)
t = 5.0
@test isapprox(postwid(s -> 1/s^4, t, v = _PWcoeffs(36)), t^3/6, rtol = 1e-12)

# test across t array
t = 1.0:10.0

setprecision(53)
@test isapprox(postwid(s -> 1/s^4, t, v = _PWcoeffs(18)), t.^3 ./6, rtol = 1e-4)

setprecision(128)
@test isapprox(postwid(s -> 1/s^4, t, v = _PWcoeffs(36)), t.^3 ./6, rtol = 1e-12)


#### F(s) = 1/(s+1), f(t) = exp(-t)
## for larger t you need a higher amount of coeffs depending on accuracy

# test for single time point
setprecision(53)
t = 1.0
@test isapprox(postwid(s -> 1/(s+1), t, v = _PWcoeffs(18)), exp(-t), rtol = 1e-4)
t = 5.0
@test isapprox(postwid(s -> 1/(s+1), t, v = _PWcoeffs(18)), exp(-t), rtol = 1e-3)

setprecision(128)
t = 1.0
@test isapprox(postwid(s -> 1/(s+1), t, v = _PWcoeffs(36)), exp(-t), rtol = 1e-12)
t = 5.0
@test isapprox(postwid(s -> 1/(s+1), t, v = _PWcoeffs(36)), exp(-t), rtol = 1e-8)

# test across t array
t = 1.0:10.0

setprecision(53)
@test isapprox(postwid(s -> 1/(s+1), t, v = _PWcoeffs(18)), exp.(-t), rtol = 1e-3)

setprecision(128)
@test isapprox(postwid(s -> 1/(s+1), t, v = _PWcoeffs(36)), exp.(-t), rtol = 1e-8)




end # module