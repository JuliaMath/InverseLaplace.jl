
function arrFcomplex(s)
    # Laplace domain
    α = complex(-0.3, 6.0)
    return [1 2;3 4] ./ (s - α)
end

function Fcomplex(s)
    # Laplace domain
    α = complex(-0.3, 6.0)
    return 1 / (s - α)
end

@testset "Array functionality" begin
    @testset "Internal functions" begin
        # Test scalar-functionality of arrcoeff, Fcomplex=(s-a)⁻¹
        let arr = try InverseLaplace._arrcoeff(Fcomplex,4,1.0,1.0);
            InverseLaplace._arrcoeff(Fcomplex,4,1.0,1.0)
            catch
                InverseLaplace._arrcoeff(Fcomplex,4,1.0,1.0)
            end, scal = InverseLaplace._wcoeff(Fcomplex,4,1.0,1.0)

                @test arr == scal
        end

        # Test for ordering of array valued coefficients (along first dimension)
        let arr = try InverseLaplace._arrcoeff(arrFcomplex,4,1.0,1.0);
            InverseLaplace._arrcoeff(arrFcomplex,4,1.0,1.0)
            catch
                InverseLaplace._arrcoeff(arrFcomplex,4,1.0,1.0)
            end, scal = InverseLaplace._wcoeff(Fcomplex,4,1.0,1.0)
                @test arr[:,1,1] == scal
                @test arr[:,1,2] == 2 .* arr[:,1,1]
                @test isapprox(arr[:,2,1] , 3 .* arr[:,1,1], atol=1E-15)
                @test arr[:,2,2] == 4 .* arr[:,1,1]
        end

        # Compare _get_coefficients and _laguerre for array functions and scalar functions
        let arr = try InverseLaplace._get_array_coefficients(arrFcomplex,4,1.0,1.0,Complex);
            InverseLaplace._get_array_coefficients(arrFcomplex,4,1.0,1.0,Complex)
            catch
                InverseLaplace._get_array_coefficients(arrFcomplex,4,1.0,1.0,Complex)
            end, scal = InverseLaplace._get_coefficients(Fcomplex,4,1.0,1.0,Complex)
                @test arr[:,1,1] == scal
                @test arr[:,1,2] == 2 .* scal
                @test isapprox(arr[:,2,1] , 3 .* scal, atol=1E-15)
                @test arr[:,2,2] == 4 .* scal

                laguerreeval = reshape(mapslices(i -> InverseLaplace._laguerre(i,1.0),arr,dims=(1)),(2,2))
                @test isapprox(laguerreeval, [1 2; 3 4] .* InverseLaplace._laguerre(scal,1.0), atol = 1E-15)
        end
    end
    @testset "ILT for arrays" begin
        # Test the ILT calculation for arrays and compare with the scalar Weeks method
        let
            # Weeks parameters
                N = 4
                σ = 1.0
                b = 0.5
                t = 1.0
                testingrank = 2
                funcdims = size(arrFcomplex(rand(Complex{Float64}))) # (2,2)

            # Evaluating ILT for array valued function
            coef = InverseLaplace._get_array_coefficients(arrFcomplex,4,σ,b,Complex)
            lag = reshape(mapslices(i -> InverseLaplace._laguerre(i,2 * b * t),coef,dims=(1)),(2,2))
            inverse = lag .* exp((σ - b) * t)

            # Evaluating ILT for scalar function, then multiplying by an array
            Ft = Weeks(Fcomplex,4,σ,b,datatype=Complex)
            Ft_array_eval = [1 2;3 4]  .* Ft(1.0)

            @test isapprox(Ft_array_eval , inverse, atol=1E-15)

            # Test the array constructor for Weeks{T}, both real and complex
            Arr_c = Weeks(arrFcomplex,N,σ,b,datatype=Complex,rank=testingrank)
            Scal_c = Weeks(Fcomplex,N,σ,b,datatype=Complex)
            Arr_r = Weeks(arrFcomplex,N,σ,b,rank=testingrank)
            Scal_r = Weeks(Fcomplex,N,σ,b)
            @test ndims(Arr_c.coefficients) == testingrank + 1
            @test size(Arr_c.coefficients) == (Arr_c.Nterms, 2,2)
            @test ndims(Arr_r.coefficients) == testingrank + 1
            @test size(Arr_r.coefficients) == (Arr_r.Nterms,2,2)
            @test real(Arr_c.coefficients) == Arr_r.coefficients
            @test isapprox(Arr_c(t),[1 2; 3 4] .* Scal_c(t),atol=1E-15)
            @test isapprox(Arr_r(t),[1 2; 3 4] .* Scal_r(t),atol=1E-15)
            @test isapprox(real(Arr_c(t)),Arr_r(t),atol=1E-15)
        end
    end
end
