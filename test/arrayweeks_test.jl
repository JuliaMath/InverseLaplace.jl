
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

# Test scalar-functionality of arrcoeff, Fcomplex=(s-a)⁻¹
let arr = try InverseLaplace._arrcoeff(Fcomplex,4,1.0,1.0);
    InverseLaplace._arrcoeff(Fcomplex,4,1.0,1.0)
    catch
        InverseLaplace._arrcoeff(Fcomplex,4,1.0,1.0)
    end, scal = InverseLaplace._wcoeff(Fcomplex,4,1.0,1.0)

        @test arr == scal
end

let arr = try InverseLaplace._arrcoeff(arrFcomplex,4,1.0,1.0);
    InverseLaplace._arrcoeff(arrFcomplex,4,1.0,1.0)
    catch
        InverseLaplace._arrcoeff(arrFcomplex,4,1.0,1.0)
    end, scal = InverseLaplace._wcoeff(Fcomplex,4,1.0,1.0)
        @test arr[:,1,1] == scal
        @test arr[:,2,1] == 2 .* arr[:,1,1]
        @test isapprox(arr[:,1,2] , 3 .* arr[:,1,1], atol=1E-15)
        @test arr[:,2,2] == 4 .* arr[:,1,1]
end
