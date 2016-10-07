using Optim

function weeks(F, t::AbstractVector, N, sig, b)
    a0 = real(wcoef(F,N,sig,b))
    a = a0[N+1:2*N]
    L = laguer(a,2*b*t)
    L .* exp((sig-b)*t)
end

function weeks(F, t, N, sig, b)
    a0 = real(wcoef(F,N,sig,b))
    a = a0[N+1:2*N]
    laguer(a,2*b*t) * exp((sig-b)*t)
end

function weekse(F, t::AbstractVector, N, sig, b)
    M = 2 * N
    a   = real(wcoef(F,M,sig,b))
    a1  = a[2*N+1:3*N]; sa1 = sum(abs(a1))
    a2  = a[3*N+1:4*N]; sa2 = sum(abs(a2))
    L   = laguer(a1,2*b*t)
    f   = L .* exp((sig-b)*t)
    est = exp(sig*t)*(sa2+eps()*sa1)
    (f,est)
end

function weekse(F, t, N, sig, b)
    M = 2 * N
    a   = real(wcoef(F,M,sig,b))
    a1  = a[2*N+1:3*N]; sa1 = sum(abs(a1))
    a2  = a[3*N+1:4*N]; sa2 = sum(abs(a2))
    L   = laguer(a1,2*b*t);
    f   = L*exp((sig-b)*t);
    est = exp(sig*t)*(sa2+eps()*sa1);
    (f,est)
end


function laguer(a::AbstractVector,x::AbstractVector)
    N = length(a) - 1
    unp1 = zeros(x)
    un = a[N+1]*ones(x)
    local unm1
    for n in N:-1:1
        unm1 =  (1//n)*(2*n-1-x) .* un - n/(n+1)*unp1 + a[n]
        unp1 = un
        un = unm1
    end
    unm1
end

function laguer(a::AbstractVector,x)
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

function wcoef(F,N,sig,b)
    n = -N:1:N-1
    h = pi/N    
    th = h*(n+1//2)
    y = b*cot(th/2)
    imunit = Complex(zero(eltype(y)),one(eltype(y)))
    s = sig +  imunit * y
    FF0 = map(F, s)
    FF = FF0 .* (b+imunit*y)
    a = fftshift(fft(fftshift(FF)))/(2*N)
    exp(Complex(zero(h),-one(h)) * n*h/2) .* a
end

function wpar2(F, t, N, sig0, sigmax, bmax)
    so = Optim.minimizer(optimize( sig -> werr2e(sig, F,t,N, sig0, sigmax, bmax), sig0, sigmax))
    bo = Optim.minimizer(optimize( b -> werr2t(b, F, N, so), 0, bmax))
    (so,bo)
end

function werr2e(sig,F,t,N,sig0,sigmax,bmax)
    b = Optim.minimizer(optimize( (b) -> werr2t(b, F,N, sig) , 0.0, bmax))
    M = 2*N
    a = wcoef(F,M,sig,b)
    a1 = a[2*N+1:3*N]
    sa1 = sum(abs(a1))
    a2 = a[3*N+1:4*N]
    sa2 = sum(abs(a2))
    sig*t + log(sa2+eps()*sa1)  # the original computed exp of this and then the log
end

function werr2t(b, F, N, sig)
    M = 2*N
    a = wcoef(F,M,sig,b)
    sa2 = sum(abs(a[3*N+1:4*N]))
    log(sa2)
end
