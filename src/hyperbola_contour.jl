### Implement Laplace transform along a hyperbola contour ###

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
    a =  zero(Complex{eltype(N)})
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
    a = zeros(Complex{eltype(h)}, Int(N))
    sk = similar(a)
    Threads.@threads for k in 0:Int(N)-1
        sk[k+1] = s_fixed((k + 1/2)*h, μ)
        dsk = ds_fixed((k + 1/2)*h, μ)
        a[k+1] = f(sk[k+1])*dsk
    end
    return a, sk, h
end

function LT_hyper_fixed(a, sk, h, t::AbstractFloat)
    b = zero(eltype(a))
    for ind in eachindex(sk)
        b += a[ind]*exp(sk[ind]*t)
    end
    return imag(b)*h/π
end

function LT_hyper_fixed(f::Function, N, t::AbstractArray)
    a, sk, h = fixed_sk(f, N, t)
    out = zeros(eltype(h), length(t))
    Threads.@threads for ind in eachindex(t)
        out[ind] = LT_hyper_fixed(a, sk, h, t[ind])
    end
    return out
end




### Implement Laplace transform along improved Talbot contour ###

## adaptive contour for each time point (can't precompute f(ω) or sk)

function z(θ, N, t; α = 0.6407, c = 1.3580)
    B = c*(sin(α*pi))^2 / ( 2*α*c^2 * (sin(α*pi))^2 - pi*sin(2*α*pi)*(sinh(α*c))^2)
    σ = 2*α*c^2 * B
    μ = 2*(sinh(α*c))^2 * B
    ν = (sinh(2*α*c) - 2*α*c)*B

    # return both z(θ) and z'(θ)
    return N/t * (-σ + μ*θ*cot(α*θ) + ν*im*θ), N/t * (ν*im  + μ*(cot(α*θ) - α*θ*(csc(α*θ))^2))
end

# compute the laplace transform along a hyperbola contour fixed time point
function talbot_improved(f::Function, N, t::AbstractFloat)
    a =  zero(Complex{eltype(N)})
    h = 2*π/N
    for k in 1:N
        zk, dzk = z(-π + (k - 1/2)*h, N, t)       
        a += f(zk)*exp(zk*t)*dzk
    end
    return a/N/im
end

function talbot_improved(f::Function, N, t::AbstractFloat)
    a =  zero(Complex{eltype(N)})
    h = 2*π/N
    for k in 1:N/2 # for real-valued f(t) need only to consider half of
        zk, dzk = z(-π + (k - 1/2)*h, N, t)       
        a += f(zk)*exp(zk*t)*dzk
        println(-π + (k - 1/2)*h)
    end
    return 2*a/N/im
end


function talbot_improved(f::Function, N, t::AbstractFloat)
    a =  zero(Complex{eltype(N)})
    h = π/N
    for k in 1:N # for real-valued f(t) need only to consider half of evaluations
        zk, dzk = z(-π + (k - 1/2)*2*h, N, t)       
        a += f(zk)*exp(zk*t)*dzk
    end
    return real(a/N/im)
end

