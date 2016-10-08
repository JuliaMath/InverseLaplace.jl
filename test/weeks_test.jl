@test typeof(Weeks(s -> 1/s)) == Weeks
@test typeof(WeeksErr(s -> 1/s)) == WeeksErr

fl = Weeks( s -> 1/s^2 )

@test_approx_eq( fl([1 2; 3 4]), [1 2; 3 4])
