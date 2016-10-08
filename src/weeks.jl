using Optim

abstract AbstractWeeks <: AbstractILt

function Base.show{T<:AbstractWeeks}(io::IO, w::T)
    print(io, string(typeof(w)), "(Nterms=", w.Nterms, ",sigma=", w.sigma, ",b=", w.b,')')
end

#### Weeks

type Weeks <: AbstractWeeks
    func::Function
    Nterms::Int
    sigma::Float64
    b::Float64
    coefficients::Array{Float64,1}
end

function _get_coefficients(func, Nterms, sigma, b)
    a0 = real(_wcoeff(func,Nterms,sigma,b))
    a = a0[Nterms+1:2*Nterms]
end

_get_coefficients(w::Weeks) = _get_coefficients(w.func, w.Nterms, w.sigma, w.b)


Weeks(func::Function, Nterms::Integer=64, sigma=1.0, b=1.0) =  Weeks(func,Nterms,sigma,b,_get_coefficients(func,Nterms,sigma,b))

function eval_ilt(w::Weeks, t)
    L = _laguerre(w.coefficients,2*w.b*t) 
    L * exp((w.sigma-w.b)*t)
end

function eval_ilt(w::Weeks, t::AbstractArray)
    L = _laguerre(w.coefficients,2*w.b*t)
    [ L1 * exp((w.sigma-w.b)*t1) for (t1,L1) in zip(t,L)]
end

function optimize(w::Weeks, t)
    (w.sigma, w.b) = _optimize_sigma_and_b(w.func, t, w.Nterms, 0.0, 30, 30)
    w.coefficients = _get_coefficients(w)
    w
end

function optimize(w::Weeks, t, N)
    w.Nterms = N
    optimize(w,t)
end

function opteval(w::Weeks, t)
    optimize(w,t)
    w(t)
end

function opteval(w::Weeks, t, N)
    optimize(w,t, N)
    w(t)
end

#### WeeksErr

type WeeksErr <: AbstractWeeks
    func::Function
    Nterms::Int
    sigma::Float64
    b::Float64
    coefficients::Array{Float64,1}
    sa1::Float64    
    sa2::Float64
end

function WeeksErr(func::Function, Nterms::Integer=64, sigma=1.0, b=1.0)
    M = 2 * Nterms
    a0 = real(_wcoeff(func,M,sigma,b))
    a1 = a0[2*Nterms+1:3*Nterms]
    sa1 = sum(abs(a1))
    sa2 = sum(abs(@view a0[3*Nterms+1:4*Nterms]))
    WeeksErr(func,Nterms,sigma,b,a1,sa1,sa2)
end

function eval_ilt(w::WeeksErr, t)
    L = _laguerre(w.coefficients,2*w.b*t) 
    f = L * exp((w.sigma-w.b)*t)
    est = exp(w.sigma*t)*(w.sa2+eps()*w.sa1)
    (f,est)
end

function eval_ilt(w::WeeksErr, t::AbstractArray)
    L = _laguerre(w.coefficients,2*w.b*t) 
    f = L .* exp((w.sigma-w.b)*t)
    est = exp(w.sigma*t)*(w.sa2+eps()*w.sa1)
    (f,est)
end

for w in (:Weeks, :WeeksErr)
    @eval (w::$(w))(t) = eval_ilt(w,t)
end

##### internal functions

function _wcoeff(F,N,sig,b)
    n = -N:1:N-1
    h = pi/N
    th = h*(n+1//2)
    y = b*cot(th/2)
    imunit = Complex(zero(eltype(y)),one(eltype(y)))
    s = sig +  imunit * y
    FF0 = map(F, s)
    FF = [ FF1 * (b + imunit * y1) for (FF1,y1) in zip(FF0,y)]
    a = fftshift(fft(fftshift(FF)))/(2*N)
    exp(Complex(zero(h),-one(h)) * n*h/2) .* a
end

function _laguerre(a::AbstractVector,x::AbstractArray)
    N = length(a) - 1
    unp1 = zeros(x)    
    un = a[N+1]*ones(x)
    local unm1
    for n in N:-1:1
#        unm1 =  (1//n)*(2*n-1-x) .* un - n/(n+1)*unp1 + a[n]
        unm1 =  [(1//n)*(2*n-1-x0) * un0 - n/(n+1)*unp10 + a[n] for (x0,un0,unp10) in zip(x,un,unp1)]
        unp1 = un
        un = unm1
    end
    unm1
end

function _laguerre(a::AbstractVector,x)
    N = length(a) - 1
    unp1 = zero(x)
    un = a[N+1]*one(x)
    local unm1
    for n in N:-1:1
        unm1 = (1//n)*(2*n-1-x) * un - n/(n+1)*unp1 + a[n]
        unp1 = un
        un = unm1
    end
    unm1
end

function _optimize_sigma_and_b(F, t, N, sig0, sigmax, bmax)
    sigma_opt = Optim.minimizer(Optim.optimize( sig -> _werr2e(sig, F,t,N, sig0, sigmax, bmax), sig0, sigmax))
    b_opt = Optim.minimizer(Optim.optimize( b -> _werr2t(b, F, N, sigma_opt), 0, bmax))
    (sigma_opt, b_opt)
end

function _werr2e(sig,F,t,N,sig0,sigmax,bmax)
    b = Optim.minimizer(Optim.optimize( (b) -> _werr2t(b, F,N, sig) , 0.0, bmax))
    M = 2*N
    a = _wcoeff(F,M,sig,b)
    a1 = @view a[2*N+1:3*N]
    sa1 = sum(abs(a1))
    a2 = @view a[3*N+1:4*N]
    sa2 = sum(abs(a2))
    sig*t + log(sa2+eps()*sa1)
end

function _werr2t(b, F, N, sig)
    M = 2*N
    a = _wcoeff(F,M,sig,b)
    sa2 = sum(abs( @view a[3*N+1:4*N]))
    log(sa2)
end

##### old function

function weeks(F, t::AbstractArray, N, sig, b)
    a0 = real(_wcoeff(F,N,sig,b))
    a = @view a0[N+1:2*N]
    L = _laguerre(a,2*b*t)
    L .* exp((sig-b)*t)
end

function weeks(F, t, N, sig, b)
    a0 = real(_wcoeff(F,N,sig,b))
    a =  @view a0[N+1:2*N]
    _laguerre(a,2*b*t) * exp((sig-b)*t)
end

function weekse(F, t::AbstractArray, N, sig, b)
    M = 2 * N
    a   = real(_wcoeff(F,M,sig,b))
    a1  = @view a[2*N+1:3*N]
    sa1 = sum(abs(a1))
    a2  = @view a[3*N+1:4*N]
    sa2 = sum(abs(a2))
    L   = _laguerre(a1,2*b*t)
    f   = L .* exp((sig-b)*t)
    est = exp(sig*t)*(sa2+eps()*sa1)
    (f,est)
end

function weekse(F, t, N, sig, b)
    M = 2 * N
    a   = real(_wcoeff(F,M,sig,b))
    a1  = @view a[2*N+1:3*N]; sa1 = sum(abs(a1))
    a2  = @view a[3*N+1:4*N]; sa2 = sum(abs(a2))
    L   = _laguerre(a1,2*b*t);
    f   = L*exp((sig-b)*t);
    est = exp(sig*t)*(sa2+eps()*sa1);
    (f,est)
end

