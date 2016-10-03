using InverseLaplace
using Base.Test

@test_approx_eq( ilt(s -> 1/s,1.0) , 1.000000000007737)
@test_approx_eq( ilt(s -> 1/s,1) , 1.0)
@test_approx_eq( gwr(s -> 1/s,1) , 1.0)
@test_approx_eq( gwr(s -> 1/s,1.0) , 1.0000000002465594)
@test_approx_eq( ilt(s -> s/(1+s^2),1) , cos(1))
@test_approx_eq( ilt(s -> 1/(1+s),1) , exp(-1))
@test_approx_eq( gwr(s -> 1/(1+s),1) , exp(-1))
@test_approx_eq( ilt(s -> 1/s^4,1//10) * 6 , 1e-3)
@test_approx_eq( gwr(s -> 1/s^4,1//10) * 6 , 0.0010000012595365085)

@test_approx_eq( ilt(s -> 1/s^3,[.1,.2]) , [0.005,0.02])
@test_approx_eq( ilt(s -> 1/s^3,[.1,.2]) , [0.005,0.02])
@test_approx_eq( map(x -> convert(Float64,x) , ilt(s -> 1/s^3,[BigFloat(1//10) , BigFloat(2//10)])), [0.005,0.02])
if Int != Int32
    @test_approx_eq( ilt(s -> 1/s^3,[1//10,2//10]) , [2882303761517117//576460752303423488, 5764607523035027//288230376151711744])
end
