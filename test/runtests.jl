using InverseLaplace
using Base.Test

@test_approx_eq( ilt(s -> 1/s,1.0) , 1.000000000007737)
@test_approx_eq( ilt(s -> 1/s,1) , 1.0)
@test_approx_eq( gwr(s -> 1/s,1) , 1.0)
@test_approx_eq( ilt(s -> s/(1+s^2),1) , cos(1))
@test_approx_eq( ilt(s -> 1/(1+s),1) , exp(-1))
@test_approx_eq( gwr(s -> 1/(1+s),1) , exp(-1))
@test_approx_eq( ilt(s -> 1/s^4,1//10) * 6 , 1e-3)
