
# Implement Laplace transfrom along a hyperbola contour

## adaptive contour for each time point (can't precompute f(ω) or sk)

function s(θ, N, t)
    μ = 4.492075287*N/t
    ϕ = 1.172104229 
    return μ + im*μ*sinh(θ + im*ϕ)
end

function ds(θ, N, t)
    μ = 4.492075287*N/t
    ϕ = 1.172104229 
    return im*μ*cosh(θ + im*ϕ)
end

# compute the laplace transform along a hyperbola contour fixed time point
function LT_hyperbola(f::Function, N, t::AbstractFloat)
    a = 0.0 + 0.0*im
    h = 1.081792140/N
    for k in 0:N-1
        sk = s((k + 1/2)*h, N, t)
        dsk = ds((k + 1/2)*h, N, t)
        a += f(sk)*exp(sk*t)*dsk
    end
    return imag(a)*h/pi
end

# loop through an entire time array with multithreads
function LT_hyperbola(f::Function, N, t::AbstractArray)
    out = similar(t)
    Threads.@threads for ind in eachindex(t)
        out[ind] = LT_hyperbola(f, N, t[ind])
    end
    return out   
end

#### fixed integration path (can precompute)
# get the fixed integration components
function LT_hyper_coef(N, t::AbstractArray; ϕ = 1.09)
    A = acosh(((π - 2*ϕ)*t[end]/t[1] + 4*ϕ - π)/((4*ϕ - π)*sin(ϕ)))
    μ = (4*π*ϕ - π^2)*N/t[end]/A
    h = A/N
    return μ, h
end

s_fixed(θ, μ; ϕ = 1.09) = μ + im*μ*sinh(θ + im*ϕ)
ds_fixed(θ, μ; ϕ = 1.09) = im*μ*cosh(θ + im*ϕ)

function fixed_sk(f::Function, N, t::AbstractArray)
    μ, h = LT_hyper_coef(N, t)
    a = Array{Complex{Float64}}(undef, N)
    sk = similar(a)
    Threads.@threads for k in 0:N-1
        sk[k+1] = s_fixed((k + 1/2)*h, μ)
        dsk = ds_fixed((k + 1/2)*h, μ)
        a[k+1] = f(sk[k+1])*dsk
    end
    return a, sk, h
end

function LT_hyper_fixed(a, sk, h, t::AbstractFloat)
    b = 0.0 + 0.0*im
    for ind in eachindex(sk)
        b += a[ind]*exp(sk[ind]*t)
    end
    return imag(b)*h/π
end


function LT_hyper_fixed(f::Function, N, t::AbstractArray)
    a, sk, h = fixed_sk(f, N, t)
    out = similar(t)
    Threads.@threads for ind in eachindex(t)
        out[ind] = LT_hyper_fixed(a, sk, h, t[ind])
    end
    return out
end