using InverseLaplace
import InverseLaplace: talbot, gwr

using Base.Test

include("weeks_test.jl")

@test_approx_eq( talbot(s -> 1/s,1.0) , 1.000000000007737)
@test_approx_eq( talbot(s -> 1/s,1) , 1.0)
@test_approx_eq( gwr(s -> 1/s,1) , 1.0)
@test_approx_eq( gwr(s -> 1/s,1.0) , 1.0000000002465594)
@test_approx_eq( talbot(s -> s/(1+s^2),1) , cos(1))
@test_approx_eq( talbot(s -> 1/(1+s),1) , exp(-1))
@test_approx_eq( gwr(s -> 1/(1+s),1) , exp(-1))

if Int != Int32 && VERSION >= v"0.5.0"
    @test_approx_eq( talbot(s -> 1/s^4,1//10) * 6 , 1e-3)
    @test_approx_eq( gwr(s -> 1/s^4,1//10) * 6 , 0.0010000012595365085)

    # At one point, talbot returned rational with rational input. Now we convert to BigFloat
    #
    # If we write our one vectorizing code the following test passes. But not if we use whatever code julia uses to do it automatically
    #    @test_approx_eq( talbot(s -> 1/s^3,[1//10,2//10]) , Rational[2882303761517117//576460752303423488, 5764607523035027//288230376151711744])

    # This test passes
    @test_approx_eq( talbot(s -> 1/s^3,[1//10,2//10]) , [BigFloat(2882303761517117//576460752303423488), BigFloat(5764607523035027//288230376151711744)])
else
    @test_approx_eq_eps( talbot(s -> 1/s^4,1//10) * 6 , 1e-3,  1e-9)
    @test_approx_eq_eps( gwr(s -> 1/s^4,1//10) * 6 , 0.0010000012595365085, 1e-9)
end

@test_approx_eq( talbot(s -> 1/s^3,[.1,.2]) , [0.005,0.02])

# gwr is often less accurate.

if VERSION >= v"0.5.0"
    @test_approx_eq_eps( gwr(s -> 1/s^3,[.1,.2]) , [0.005,0.02], 1e-10)
else
    @test_approx_eq_eps( gwr(s -> 1/s^3,[.1,.2]) , [0.005,0.02], 1e-7)
end
   
@test_approx_eq( InverseLaplace.talbotarr(s -> 1/s^3,[.1,.2]) , [0.005,0.02])

@test_approx_eq( map(x -> convert(Float64,x) , talbot(s -> 1/s^3,[BigFloat(1//10) , BigFloat(2//10)])), [0.005,0.02])

@test typeof(talbot(s -> 1/s^3,1//10)) == BigFloat
@test typeof(talbot(s -> 1/s^3,1//10, 64)) == BigFloat
