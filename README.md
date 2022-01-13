# GeneralizedSylvesterSolver.jl

[WORK IN PROGRESS] This Julia package solves
```
a x + b x (c ⊗ c ⊗ ... ⊗ c) = d
```
by using
```
(I + c^T ⊗ c^T ⊗ ... ⊗ c^T ⊗ b)x = d
```

## Installation
The package requires Julia 1.6.3 or higher
```
using Pkg
Pkg.add("GeneralizedSylvesterSolver")
```

## Usage

```
ws = IPlusAtKronBWs(ma, mb, mc, order)
generalized_sylvester_solver!(a::AbstractMatrix, b::AbstractMatrix, c::AbstractMatrix,
                              d::AbstractMatrix, order::Int64, ws::IPlusAtKronBWs)
```
whith
 - `a` is a ma x na matrix
 - `b` is a mb x nb matrix
 - `c` is a mc x nc matrix
 - `d` is a md x nd matrix
 - `order` is an integer representing the number of occurences of `c^T` in the Kronecker products
 - `ws` is an instance of the `IPlusAtKronBWs` type   

## Version
  - 0.1.1
    
## References
O. Kamenik (2005), "Solving SDGE models: A new algorithm for the Sylvester
  equation", <i>Computational Economics 25</i>, 167--187.
