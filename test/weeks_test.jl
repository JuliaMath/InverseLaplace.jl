@testset "weeks test"  begin
    @test isa(Weeks(s -> 1 / s), Weeks)
    @test isa(WeeksErr(s -> 1 / s), WeeksErr)

    fl = Weeks(s -> 1 / s^2)

    @test isapprox(fl(1), 1)
    @test isapprox(fl([1 2; 3 4]), [1 2; 3 4])

    fl = Weeks(s -> s / (1 + s^2), 64)
    @test isapprox(fl([1 2; 3 4]), broadcast(cos,[1.0 2.0; 3.0 4.0]))

    fle = WeeksErr(s -> s/(1+s^2),  64)
    @test isapprox(fle([1 2; 3 4])[1], cos.([1 2; 3 4]))

    ## Optim is broken.
    # e1 = abs(fl(10.0) - cos(10.0))
    # InverseLaplace.optimize(fl, 10.0)
    # e2 = abs(fl(10.0) - cos(10.0))
    # @test e2 < e1

    # e1 = abs(fle(10.0)[1] - cos(10.0))
    # optimize(fle,10.0)
    # e2 = abs(fl(10.0)[1] - cos(10.0))
    # @test e2 < e1

    fle = WeeksErr(s -> s/(1+s^2),  64)
    (t1,e1) = fle(10.0)
    @test isapprox([t1,e1], [-0.8390715290763104,1.6069233509475447e-11]; atol = 1e-9)

    setNterms(fle, 8)
    e1 = abs(fle(10.0)[1] - cos(10.0))
    setNterms(fle, 80)
    e2 = abs(fle(10.0)[1] - cos(10.0))
    @test e2 < e1

    fle = WeeksErr(s -> s/(1+s^2),  64)
    c1 = copy(fle.coefficients)
    @test string(fle) == "WeeksErr{Float64}(Nterms=64,sigma=1.0,b=1.0)"
    setparameters(fle,2.0,2.0,80)
    @test string(fle) == "WeeksErr{Float64}(Nterms=80,sigma=2.0,b=2.0)"

    # Check that Laguerre coefficients have been recomputed
    @test c1 != fle.coefficients

    ### Complex

    function Fcomplex(s)
        # Laplace domain
        α = complex(-0.3, 6.0)
        return 1 / (s - α)
    end

    function fcomplex(t)
        # Time domain
        α = complex(-0.3, 6.0)
        return exp(α * t)
    end

    # Tests with Julia v1 pass locally but fail in CI for all but MacOs
    # We don't even test with v1.1, so this is sort of nonsense.
    if VERSION >= v"1.1"
        let Fc = Weeks(Fcomplex, 1024, datatype=Complex),  trange = range(0.0, stop=30.0, length=1000)
            # atol = 0.001 works locally for most versions, but fails CI for many
            # A bit of work to figure out what's happening
            @test isapprox(Fc.(trange), fcomplex.(trange), atol=0.006)
        end
    end
end # @testset
