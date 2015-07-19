# ILT

Inverse Laplace transform

```julia
    ilt(s -> 1/s^2, 0.5)
```

Use more terms for BigFloat input

```julia
    ilt(s -> 1/s^2, BigFloat(1//2), 64)
```
	
