fl = ILt(s -> 1/s^3, talbot)
@test_approx_eq_eps( fl(1.0), 0.5, 1e-9)
@test_approx_eq_eps( fl([1.0,2.0]), [0.5, 2.0] , 1e-9)

fl = ILt(s -> 1/s^3, gwr,8)
@test_approx_eq_eps( fl(1.0), 0.5, 1e-5)

p = ILtPair(Weeks(s -> s/(1+s^2)), cos)

e1 = abserr(p, 10.0)
optimize(p,10.0)
e2 = abserr(p, 10.0)
@test e2 < e1

@test typeof(Weeks(iltpair_power(5))) <: ILtPair
