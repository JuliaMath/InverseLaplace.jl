@test typeof(Weeks(s -> 1/s)) == Weeks
@test typeof(WeeksErr(s -> 1/s)) == WeeksErr

fl = Weeks( s -> 1/s^2 )

@test_approx_eq( fl(1), 1)
@test_approx_eq( fl([1 2; 3 4]), [1 2; 3 4])

fl = Weeks( s -> s/(1+s^2),  64 )
@test_approx_eq( fl([1 2; 3 4]), cos([1 2; 3 4]))

fle = WeeksErr( s -> s/(1+s^2),  64 )
@test_approx_eq( fle([1 2; 3 4])[1], cos([1 2; 3 4]))

e1 = abs(fl(10.0) - cos(10.0))
optimize(fl,10.0)
e2 = abs(fl(10.0) - cos(10.0))
@test e2 < e1

e1 = abs(fle(10.0)[1] - cos(10.0))
optimize(fle,10.0)
e2 = abs(fl(10.0)[1] - cos(10.0))
@test e2 < e1

fle = WeeksErr( s -> s/(1+s^2),  64 )
(t1,e1) = fle(10.0)
@test_approx_eq_eps([t1,e1], [-0.8390715290763104,1.6069233509475447e-11], 1e-9)

setNterms(fle, 8)
e1 = abs(fle(10.0)[1] - cos(10.0))
setNterms(fle, 80)
e2 = abs(fle(10.0)[1] - cos(10.0))
@test e2 < e1

