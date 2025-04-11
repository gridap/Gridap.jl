"""
This module provides a collection of uni- and multi-variate scalar- and multi-value'd polynomial bases.

Most of the basis polynomials are composed using products of 1D polynomials,
represented by the type [`Polynomial`](@ref).
Five `Polynomial` families are currently implemented: [`Monomial`](@ref),
[`Legendre`](@ref), [`Chebyshev`](@ref), [`Bernstein`](@ref) and [`ModalC0`](@ref).

The polynomial bases all subtype [`PolynomialBasis`](@ref), which subtypes
`AbstractVector{<:Field}`, so they implement the `Field` interface up to first
or second derivatives.

Constructors for commonly used bases (see the documentation for the spaces definitions):
- ℚ spaces: `[Polynomial]Basis(Val(D), V, order)`
- ℙ spaces: `[Polynomial]Basis(..., Polynomials._p_filter)`
- ℚₙ\\ℚₙ₋₁: `[Polynomial]Basis(..., Polynomials._qh_filter)`
- ℙₙ\\ℙₙ₋₁: `[Polynomial]Basis(..., Polynomials._ph_filter)`
- ℕ𝔻(△): [`PGradBasis`](@ref)`(Val(D), T, order)`
- ℕ𝔻(□): [`QGradBasis`](@ref)`(...)`
- ℝ𝕋(△): [`PCurlGradBasis`](@ref)`(...)`
- ℝ𝕋(□): [`QCurlGradBasis`](@ref)`(...)`

### Examples

```julia
using Gridap
using Gridap.Polynomials
using Gridap.Fields: return_type

# Basis of ℚ¹₂ of Float64 value type based on Bernstein polynomials:
#   {(1-x)², 2x(1-x), x²}
D = 1; n = 2 # spatial dimension and order
b = BernsteinBasis(Val(D), Float64, n)

# APIs
length(b)        # 3
return_type(b)   # Float64
get_order(b)     # 2

xi =Point(0.1)
evaluate(b, xi)
evaluate(Broadcasting(∇)(b), xi)  # gradients
evaluate(Broadcasting(∇∇)(b), xi) # hessians, not all basis support hessians
evaluate(b, [xi, xi]) # evaluation on arrays of points

# Basis of ℚ²₂ of Float64 value type based on Legendre polynomials, our 1D
# Legendre polynomials are normalized for L2 scalar product and moved from
# [-1,1] to [0,1] using the change of variable x -> 2x-1
#   { 1,            √3(2x-1),             √5(6x²-6x+2),
#     √3(2y-1),     √3(2x-1)√3(2y-1),     √5(6x²-6x+2)√3(2y-1),
#     √5(6y²-6y+2), √3(2x-1)√5(6x²-6x+2), √5(6x²-6x+2)√5(6y²-6y+2) }
D = 2; n = 2
b = LegendreBasis(Val(D), Float64, n)

# Basis of (ℙ³₁)³ of VectorValue{3,Float64} value type, based on monomials:
#   {(1,0,0), (0,1,0), (0,0,1)
#    (x,0,0), (0,x,0), (0,0,x)
#    (y,0,0), (0,y,0), (0,0,y)
#    (z,0,0), (0,z,0), (0,0,z)}
D = 3; n = 1
b = MonomialBasis(Val(D), VectorValue{D,Float64}, n, Polynomials._p_filter)
evaluate(b, Point(.1, .2, .3)

# a basis for Nedelec on tetrahedra with curl in ℙ₂
b = PGradBasis(Monomial, Val(3), Float64, 2)          # basis of order 3

# a basis for Nedelec on hexahedra with divergence in ℚ₂
b = QGradBasis(Bernstein, Val(3), Float64, 2)         # basis of order 3

# a basis for Raviart-Thomas on tetrahedra with divergence in ℙ₂
b = PCurlGradBasis(Chebyshev{:T}, Val(3), Float64, 2) # basis of order 3

# a basis for Raviart-Thomas on rectangles with divergence in ℚ₃
b = QCurlGradBasis(Bernstein, Val(2), Float64, 3)     # basis of order 4
```

$(public_names_in_md(@__MODULE__))
"""
module Polynomials

using DocStringExtensions
using LinearAlgebra: mul!
using LinearAlgebra: I
using StaticArrays
using Gridap.Helpers
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields

using PolynomialBases: jacobi, jacobi_and_derivative
using Combinatorics: multiexponents, multinomial, combinations
using Base.Iterators: take

import Gridap.Fields: evaluate!
import Gridap.Fields: return_cache
import Gridap.Arrays: return_type
import Base.@propagate_inbounds

export Polynomial
export isHierarchical
export Monomial
export Legendre
export Chebyshev
export ModalC0
export Bernstein

export PolynomialBasis
export get_order

export CartProdPolyBasis
export get_exponents
export get_orders
export MonomialBasis
export LegendreBasis
export ChebyshevBasis
export BernsteinBasis
export BernsteinBasisOnSimplex

export CompWiseTensorPolyBasis
export QGradBasis
export QCurlGradBasis

export NedelecPolyBasisOnSimplex
export PGradBasis

export RaviartThomasPolyBasis
export PCurlGradBasis

export ModalC0Basis

# deprecated
export num_terms
export QGradMonomialBasis
export QCurlGradMonomialBasis
export PCurlGradMonomialBasis


include("PolynomialInterfaces.jl")

include("CartProdPolyBases.jl")

include("CompWiseTensorPolyBases.jl")

include("NedelecPolyBases.jl")

include("RaviartThomasPolyBases.jl")

include("MonomialBases.jl")

include("LegendreBases.jl")

include("ChebyshevBases.jl")

include("BernsteinBases.jl")

include("ModalC0Bases.jl")

include("Combinations.jl")

include("PLambdaBases.jl")

include("Deprecated.jl")

end # module
