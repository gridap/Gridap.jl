module BezierRefFEsTests

using Gridap

using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.Arrays
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.CellData
using FillArrays
using Test


struct Bezier <: ReferenceFEName end

const bezier = Bezier()

struct BezierRefFE{D} <: ReferenceFE{D}
  reffe::ReferenceFe{D}
  node_permutations::Vector{Int32}
  shapefuns::AbstractVector{<:Field}
  metadata
end

function BezierRefFE(T::Type{T},p::Polytope{D},orders) where {D,T}
  reffe = LagrangianRefFE(T,p,orders)
  node_permutations = bezier_node_permutations(p,orders)
  shapefuns = berstein_basis(reffe.reffe.prebasis)
  lagreffaces = reffe.reffe.metadata
  reffaces = [ BezierRefFE(T,get_polytope(rf),get_orders) for rf in lagreffaces ]
  tuple( reffaces ... )
  BezierRefFE(reffe,node_permutations,shapefuns,reffaces)
end

get_node_permutations(reffe::BezierRefFE) = reffe.node_permutations

get_metadata(reffe::BezierRefFE) = reffe.metadata

get_shapefuns(reffe::BezierRefFE) = reffe.shapefuns

function compute_shapefuns(reffe::BezierRefFE)
  berstein_basis(prebasis)
end

function compute_shapefuns(reffe::BezierRefFE,weights)
  rational_berstein_basis(prebasis,weights)
end

function _bernstein_term(p,a,i,::Val{1})
  @assert i ≤ p
  @assert a ≤ p
  if ( i ≤ a ≤ p )
    f = factorial(p) ÷ ( factorial(i)*factorial(p-i) )
    b1 = binomial( p-i, a-i )
    s = (-1)^(a-i)
    f*b1*s
  else
    0
  end
end

function _bernstein_term(p,a,b,i,j,::Val{2})
  @assert i+j ≤ p
  @assert a+b ≤ p
  if ( i ≤ a ≤ p-j ) && ( j ≤ b ≤ p-a )
    f = factorial(p) ÷ ( factorial(i)*factorial(j)*factorial(p-i-j) )
    b1 = binomial( p-i-j, a-i )
    b2 = binomial( p-a-j, b-j )
    s = (-1)^(a+b-i-j)
    f*b1*b2*s
  else
    0
  end
end

function _bernstein_term(p,a,b,c,i,j,k,::Val{3})
  @assert i+j+k ≤ p
  @assert a+b+c ≤ p
  if ( i ≤ a ≤ p-j-k ) && ( j ≤ b ≤ p-a-k ) && ( k ≤ c ≤ p-a-b )
    p!,i!,j!,k! = factorial(p),factorial(i),factorial(j),factorial(k)
    f = p! ÷ ( i!*j!*k!*factorial(p-i-j-k) )
    b1 = binomial( p-i-j-k, a-i )
    b2 = binomial( p-a-j-k, b-j )
    b3 = binomial( p-a-b-k, c-k )
    s = (-1)^(a+b+c-i-j-k)
    f*b1*b2*b3*s
  else
    0
  end
end

function _bernstein_term(p,a::NTuple{<:N},i::NTuple{<:N}) where N
  _bernstein_term(p,a...,i...,Val{N}())
end

function _berstein_matrix(prebasis)
  p = get_order(prebasis)
  e = get_exponents(prebasis)
  M = CartesianIndices( (length(e),length(e)) )
  lazy_map( i-> _bernstein_term(p,e[i.I[1]],e[i.I[2]]), M)
end

function berstein_basis(prebasis)
  C = _berstein_matrix(prebasis)
  linear_combination( C, prebasis )
end

function berstein_basis_hex(prebasis::MonomialBasis{2,T}) where {T}
  o = get_orders(prebasis)
  e = get_exponents(prebasis)
  #TODO: avoid to create new 1D basis each time
  p1 = MonomialBasis{1}(T,o[1])
  p2 = MonomialBasis{1}(T,o[2])
  b1 = berstein_basis(p1)
  b2 = berstein_basis(p2)
  lazy_map( i->b1[i[1]]*b2[i[2]], e)
end

function berstein_basis_hex(prebasis::MonomialBasis{3,T}) where {T}
  o = get_orders(prebasis)
  e = get_exponents(prebasis)
  p1 = MonomialBasis{1}(T,o[1])
  p2 = MonomialBasis{1}(T,o[2])
  p3 = MonomialBasis{1}(T,o[3])
  b1 = berstein_basis(p1)
  b2 = berstein_basis(p2)
  b3 = berstein_basis(p3)
  lazy_map( i->b1[i[1]]*b2[i[2]]*b3[i[3]], e)
end

function rational_bernstein_basis(prebasis,weights)
  #TODO: `/` operator not working for Type{Fields}
  basis = berstein_basis(prebasis)
  basis ./ ( weights ⋅ basis )
end

function bezier_node_permutations(p::Polynomials,orders)
  @notimplemented
  #TODO: general or hard-coded?
  # use _coord_to_terms(coords,orders), or similar
end

p_filter(e,o) = sum(e) ≤ o



# CURVE
p = 2

prebasis_seg =MonomialBasis{1}(Float64,p,p_filter)

C = _berstein_matrix(prebasis_seg)
C12 =
[
 1  -2   1
 0   2  -2
 0   0   1
]

@test transpose(C) == C12

ϕ = berstein_basis(prebasis_seg)

X = [ Point(0.0,0.0), Point(0.5,0.5), Point(1.0,0.0) ]

ψ = linear_combination(X,ϕ)

ξ = [Point(0.0),Point(0.5),Point(1.0)]

Ψ = Fill(ψ,length(ξ))

Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == [ Point(0.0,0.0), Point(0.5,0.25), Point(1.0,0.0) ]


p = 3
prebasis_seg =MonomialBasis{1}(Float64,p,p_filter)

C = _berstein_matrix(prebasis_seg)

C13 =
[
 1  -3   3  -1
 0   3  -6   3
 0   0   3  -3
 0   0   0   1
]

@test transpose(C) == C13

ϕ = berstein_basis(prebasis_seg)

X = [ Point(0.0,0.0), Point(0.25,0.5), Point(0.75,0.5), Point(1.0,0.0) ]

ψ = linear_combination(X,ϕ)

ξ = [Point(0.0),Point(0.5),Point(1.0)]

Ψ = Fill(ψ,length(ξ))

Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == [ Point(0.0,0.0), Point(0.5,0.375), Point(1.0,0.0) ]

# TRIANGLE

p = 2

prebasis_tri =MonomialBasis{2}(Float64,p,p_filter)

C = _berstein_matrix(prebasis_tri)

C22 =
[
 1  -2   1  -2   2   1
 0   2  -2   0  -2   0
 0   0   1   0   0   0
 0   0   0   2  -2  -2
 0   0   0   0   2   0
 0   0   0   0   0   1
]

@test transpose(C) == C22

ϕ = berstein_basis(prebasis_tri)

X = [ Point(0.0,0.0), Point(0.5,0.0), Point(1.0,0.0), Point(0.0,0.5), Point(0.5,0.5), Point(0.0,1.0) ]

ψ = linear_combination(X,ϕ)

ξ = [Point(0.0,0.0),Point(0.0,0.5),Point(0.3,0.3),Point(1.0,0.0),Point(0.0,1.0)]

Ψ = Fill(ψ,length(ξ))

Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == ξ

p = 3

prebasis_tri =MonomialBasis{2}(Float64,p,p_filter)

C = _berstein_matrix(prebasis_tri)

C23 =
[
 1  -3   3  -1  -3   6  -3   3  -3  -1
 0   3  -6   3   0  -6   6   0   3   0
 0   0   3  -3   0   0  -3   0   0   0
 0   0   0   1   0   0   0   0   0   0
 0   0   0   0   3  -6   3  -6   6   3
 0   0   0   0   0   6  -6   0  -6   0
 0   0   0   0   0   0   3   0   0   0
 0   0   0   0   0   0   0   3  -3  -3
 0   0   0   0   0   0   0   0   3   0
 0   0   0   0   0   0   0   0   0   1
]

@test transpose(C) == C23


X = [
  Point(0.0,0.0), Point(0.25,0.0), Point(0.75,0.0), Point(1.0,0.0),
  Point(0.0,0.25), Point(0.25,0.25), Point(0.75,0.25),
  Point(0.0,0.75), Point(0.25,0.75),
  Point(0.0,1.0) ]

ϕ = berstein_basis(prebasis_tri)

ψ = linear_combination(X,ϕ)

ξ = [Point(0.0,0.0),Point(0.0,0.5),Point(0.5,0.5),Point(1.0,0.0),Point(0.0,1.0)]

Ψ = Fill(ψ,length(ξ))

Xi = lazy_map( evaluate, Ψ, ξ )

@test Xi == ξ


#TODO: Ready to build reffe
#     Integrate w/ quad points as Tutorial13


end # module
