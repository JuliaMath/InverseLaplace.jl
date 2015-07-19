module ILT

export ilt

# Inverse Laplace transform
# This seems to work correctly for big float input.
function ilt(func, t, M)
    bM = convert(typeof(t),M)    
    r = (2 * bM) / (5*t)
    term = (1//2) * exp(r*t) * func(r)
    for i in 1:M-1
        theta = i * (pi/bM)
        s = r*theta*(complex(cot(theta),one(theta)))
        sigma = theta + (theta*cot(theta)-1)*cot(theta)
        term += real(exp(t*s) * complex(one(t),sigma) * func(s))
    end
    return term * 2 / (5*t)
end

# 32 gives more stable results than 64 for some reason.
# Well, depends on precision of floating point ops.
# 32 works better for machine precision
function ilt(func,t)
    ilt(func,t,32)
end

end # module
