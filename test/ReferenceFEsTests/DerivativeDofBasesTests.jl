module DofsTests

using Test
using Gridap
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.Arrays

using LinearAlgebra: I
using Gridap.ReferenceFEs: PointDerivativeValue, PointDerivativeDofBasis
import Gridap.Arrays: evaluate, return_cache, evaluate!

function return_cache(σ::PointDerivativeValue{P}, f) where P
  println("caching $σ $f")
end

function evaluate!(c, σ::PointDerivativeValue{P}, f) where P
  x = σ.point
  ∇f = Broadcasting(∇)(f)
  ∇fx = evaluate(∇f,x)
  σf = ∇fx⋅σ.direction
end

function evaluate!(c, σ::PointDerivativeValue{P}, f::AbstractArray{<:Field}) where P
  x = σ.point
  u = σ.direction
  ∇f = Broadcasting(∇)(f)
  ∇f_dot_u = Broadcasting(Operation(⋅))(∇f, u) # We have to extend the basis to 3D
  σf = evaluate(∇f_dot_u, x)
  #σf = ∇fx⋅σ.direction
end


D, T = 2, Float64
b = FEEC_poly_basis(Val(D),T,3,0,:P⁻,Monomial)
f = linear_combination( collect( i==2 for i in eachindex(b)), b)

P = Point{D,T}
σ  = PointDerivativeValue{P}(P(0,0), P(1,0))
σ2 = PointDerivativeValue{P}(P(0,0), P(0,1))

evaluate(f, σ.point)
∇f = ∇(f)
evaluate(∇f, σ.point)

c = return_cache(σ, f)
evaluate!(c, σ, f)
evaluate!(c, σ2, f)

f = linear_combination(I(length(b)), b)
∇f = Broadcasting(∇)(f)
evaluate(∇f, σ.point)

∇f = Broadcasting(Operation(∇))(f)
evaluate(∇f, σ.point)

c = return_cache(σ, f)
evaluate!(c, σ, f)
evaluate!(c, σ2, f)

prebasis = FEEC_poly_basis(Val(1), T, 3, 0,:P, Monomial) # PᵣΛᴰ⁻¹, r = order
ReferenceFEs._hermite_dofs(T, SEGMENT, prebasis)

end # module
