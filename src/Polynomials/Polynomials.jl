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
- ℚ  spaces: `[Polynomial]Basis(Val(D), V, order)`
- ℙ  spaces: `[Polynomial]Basis(..., Polynomials._p_filter)`
- 𝕊r spaces: `[Polynomial]Basis(..., Polynomials._ser_filter)`
- ℚ̃  spaces: `[Polynomial]Basis(..., Polynomials._qh_filter)`
- ℙ̃  spaces: `[Polynomial]Basis(..., Polynomials._ph_filter)`

For bases for the Nélélec, Raviart-Thomas and BDM element spaces, use
[`FEEC_poly_basis`](@ref) with the arguments found in the
[ReferenceFEs summary](@ref "Reference FE summary") of the documentation.

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

# a basis for Nedelec on tetrahedra with curl in ℙ³₂
D, k, r = 3, 1, 2+1
b = FEEC_poly_basis(Val(D),Float64,r,k,:P⁻)                 # basis of order 3

# a basis for Nedelec on hexahedra with curl in ℚ³₁
D, k, r = 3, 1, 1+1
b = FEEC_poly_basis(Val(D),Float64,r,k,:Q⁻)                 # basis of order 2

# a basis for Raviart-Thomas on quadrilateral with divergence in ℚ₁
D, k, r = 2, 2-1, 1+1
b = FEEC_poly_basis(Val(D),Float64,r,k,:Q⁻; rotate_90=true) # basis of order 3

# a basis for Raviart-Thomas on tetrahedra with divergence in ℙ₂
D, k, r = 3, 3-1, 2+1
b = FEEC_poly_basis(Val(D),Float64,r,k,:P⁻)                 # basis of order 3
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
using Base: @propagate_inbounds

import Gridap.Fields: evaluate!
import Gridap.Fields: return_cache
import Gridap.Arrays: return_type
import Gridap.Arrays: testvalue

export Polynomial
export isHierarchical

export PolynomialBasis
export get_order

export CartProdPolyBasis
export get_exponents
export get_orders
export Monomial
export MonomialBasis
export LegendreBasis
export Legendre
export ChebyshevBasis
export Chebyshev
export Bernstein
export BernsteinBasis

export BernsteinBasisOnSimplex
export bernstein_terms
export bernstein_term_id

export CompWiseTensorPolyBasis
export NedelecPolyBasisOnSimplex
export RaviartThomasPolyBasis

export ModalC0Basis
export ModalC0

export PLambdaBasis
export PmLambdaBasis
export PmΛ_bubbles
export PΛ_bubbles
export get_bubbles
export print_indices
#export get_FEEC_poly_degree
#export get_FEEC_form_degree
#export get_FEEC_family

export FEEC_space_definition_checks
export FEEC_poly_basis


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

include("PLambdaBases.jl")

include("ExteriorCalculusBases.jl")

include("Deprecated.jl")

end # module
