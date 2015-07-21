# ILT

Numerical inverse Laplace transform

Example. Evaluate the inverse Laplace transform of `1/s^2` at
`t=0.5`.
```julia
    ilt(s -> 1/s^2, 0.5)
```

The third argument gives the number of terms used.
Use more terms for BigFloat input
```julia
    ilt(s -> 1/s^2, BigFloat(1//2), 64)
```
