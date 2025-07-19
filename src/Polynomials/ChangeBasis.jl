
"""
    change_basis(basis,changeofbasis::AbstractMatrix)

# Examples

Compute the Lagrangian basis associated with a set of nodes

```jldoctests
using LinearAlgebra
using Gridap.Fields
using Gridap.Polynomials

D = 2
order = 1
f = MonomialBasis(Val(D),Float64,order)

nodes = Point{2,Int}[(0,0),(1,0),(0,1),(1,1)]
change = inv(evaluate(f,nodes))

g = change_basis(f,change)
println(evaluate(g,nodes))

# output
[1.0 0.0 0.0 0.0; 0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0]

```
"""
function change_basis(basis,changeofbasis::AbstractMatrix)
  BasisFromChangeOfBasis(basis,changeofbasis)
end

function return_type(::typeof(change_basis),prebasis,matrix_inv)
  typeof(change_basis(prebasis,matrix_inv))
end

struct BasisFromChangeOfBasis{B,M} <: AbstractVector{Field}
  basis::B
  change::M
  function BasisFromChangeOfBasis(basis,change::AbstractMatrix)
    B = typeof(basis)
    M = typeof(change)
    new{B,M}(basis,change)
  end
end

struct BasisTermFromChangeOfBasis end

Base.size(a::BasisFromChangeOfBasis) = (length(a.basis),)
Base.axes(a::BasisFromChangeOfBasis) = (axes(a.basis,1),)
# @santiagobadia : Not sure we want to create the real computation here
Base.getindex(a::BasisFromChangeOfBasis,i::Integer) = BasisTermFromChangeOfBasis()
Base.IndexStyle(::BasisFromChangeOfBasis) = IndexLinear()

function return_cache(b::BasisFromChangeOfBasis,x)
  cb = return_cache(b.basis,x)
  bx = evaluate!(cb,b.basis,x)
  c = CachedArray(bx*b.change)
  (c,cb)
end

function evaluate!(cache,b::BasisFromChangeOfBasis,x)
  c, cb = cache
  bx = evaluate!(cb,b.basis,x)
  setsize!(c,size(bx))
  mul!(c.array,bx,b.change)
  c.array
end

# Aren't next 4 functions out of date and removable ?
function return_gradient_cache(b::BasisFromChangeOfBasis,x)
  cb = return_gradient_cache(b.basis,x)
  bx = evaluate_gradient!(cb,b.basis,x)
  c = CachedArray(bx*b.change)
  (c,cb)
end

function evaluate_gradient!(cache,b::BasisFromChangeOfBasis,x)
  c, cb = cache
  bx = evaluate_gradient!(cb,b.basis,x)
  setsize!(c,size(bx))
  mul!(c.array,bx,b.change)
  c.array
end

function return_hessian_cache(b::BasisFromChangeOfBasis,x)
  cb = return_hessian_cache(b.basis,x)
  bx = evaluate_hessian!(cb,b.basis,x)
  c = CachedArray(bx*b.change)
  (c,cb)
end

function evaluate_hessian!(cache,b::BasisFromChangeOfBasis,x)
  c, cb = cache
  bx = evaluate_hessian!(cb,b.basis,x)
  setsize!(c,size(bx))
  mul!(c.array,bx,b.change)
  c.array
end
