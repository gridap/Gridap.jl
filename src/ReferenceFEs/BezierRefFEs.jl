"""
    struct Bezier  <: ReferenceFEName
"""
struct Bezier <: ReferenceFEName end

"""
    const bezier = Bezier()

Singleton of the [`Bezier`](@ref) reference FE name.
"""
const bezier = Bezier()

struct BezierRefFE{D} <: LagrangianRefFE{D}
  reffe::ReferenceFE{D}
  node_to_own_node::Vector{Int}
  shapefuns::AbstractVector{<:Field}
  metadata
end

function BezierRefFE(::Type{T},p::Polytope{D},orders) where {D,T}
  reffe = LagrangianRefFE(T,p,orders)
  nodes = get_node_coordinates(reffe)
  prebasis = get_prebasis(reffe)
  node_to_own_node = compute_node_to_bezier_node(prebasis,nodes)
  shapefuns = berstein_basis(prebasis,p)
  shapefuns = collect(lazy_map(Reindex(shapefuns),node_to_own_node))
  lagreffaces = reffe.reffe.metadata
  reffaces = _convert_reffaces(T,lagreffaces)
  tuple( reffaces ... )
  BezierRefFE(reffe,node_to_own_node,shapefuns,reffaces)
end

function ReferenceFE(
  polytope::Polytope,
  ::Bezier,
  ::Type{T},
  orders::Union{Integer,Tuple{Vararg{Integer}}}) where T

  BezierRefFE(T,polytope,orders)
end

get_node_to_own_node(reffe::BezierRefFE) = reffe.node_to_own_node

get_metadata(reffe::BezierRefFE) = reffe.metadata

get_shapefuns(reffe::BezierRefFE) = reffe.shapefuns

get_prebasis(reffe::BezierRefFE) = get_prebasis(reffe.reffe)

get_node_coordinates(reffe::BezierRefFE) = get_node_coordinates(reffe.reffe)

get_vertex_node(reffe::BezierRefFE) = get_vertex_node(reffe.reffe)

get_face_nodes(reffe::BezierRefFE) = get_face_nodes(reffe.reffe)

get_polytope(reffe::BezierRefFE) = get_polytope(reffe.reffe)

get_orders(reffe::BezierRefFE) = get_orders(reffe.reffe)

get_order(reffe::BezierRefFE) = get_order(reffe.reffe)

is_first_order(reffe::BezierRefFE) = is_first_order(reffe.reffe)

Conformity(reffe::BezierRefFE) = Conformity(reffe.reffe)

get_dof_basis(reffe::BezierRefFE) = get_dof_basis(reffe.reffe)

function get_face_own_nodes(reffe::BezierRefFE,conf::GradConformity)
  get_face_own_nodes(reffe.reffe,conf)
end

get_face_dofs(reffe::BezierRefFE) = get_face_dofs(reffe.reffe)

num_dofs(reffe::BezierRefFE) = num_dofs(reffe.reffe)

function (==)(a::BezierRefFE{D},b::BezierRefFE{D}) where D
  t = true
  a.reffe == b.reffe
  na = get_node_to_own_node(a)
  nb = get_node_to_own_node(a)
  t = t && (na == nb)
  t
end

function (==)(a::BezierRefFE,b::BezierRefFE)
  false
end

function _convert_reffaces(::Type{T},reffaces) where T
  _reffaces = []
  for d_reffaces in reffaces
    _d_reffaces = [ _lagrange_to_bezier_fe(T,reffe) for reffe in d_reffaces ]
    push!(_reffaces,_d_reffaces)
  end
  tuple(_reffaces...)
end

function _lagrange_to_bezier_fe(::Type{T},reffe::LagrangianRefFE) where T
  BezierRefFE(T,get_polytope(reffe),get_orders(reffe))
end

function compute_node_to_bezier_node(prebasis::MonomialBasis{D,T},nodes) where {D,T}
  orders = get_orders(prebasis)
  terms = _coords_to_terms(nodes,orders)
  _prebasis = MonomialBasis(Val(D),T,orders,terms)
  _exps = get_exponents(_prebasis)
  exps = get_exponents(prebasis)
  [ findfirst( isequal(i), exps) for i in _exps ]
end

## Bernstein Basis

function _bernstein_term(p,a,i)
  @check i ≤ p
  @check a ≤ p
  if ( i ≤ a ≤ p )
    f = factorial(p) ÷ ( factorial(i)*factorial(p-i) )
    b1 = binomial( p-i, a-i )
    s = (-1)^(a-i)
    f*b1*s
  else
    0
  end
end

function _bernstein_term(p,a,b,i,j)
  @check i+j ≤ p
  @check a+b ≤ p
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

function _bernstein_term(p,a,b,c,i,j,k)
  @check i+j+k ≤ p
  @check a+b+c ≤ p
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

function _bernstein_term_n_cube(p,a,b,i,j)
  _bernstein_term(p,a,i)*_bernstein_term(p,b,j)
end

function _bernstein_term_n_cube(p,a,b,c,i,j,k)
  _bernstein_term(p,a,i)*_bernstein_term(p,b,j)*_bernstein_term(p,c,k)
end

function _bernstein_term(
  p::Polytope,
  order::Integer,
  a::NTuple{N},
  i::NTuple{N}) where N

  args = (order,a...,i...)
  if is_simplex(p)
    _bernstein_term(args...)
  else
    _bernstein_term_n_cube(args...)
  end
end

function _berstein_matrix(prebasis,p::Polytope)
  o = get_order(prebasis)
  e = get_exponents(prebasis)
  M = CartesianIndices( (length(e),length(e)) )
  collect(lazy_map( i-> _bernstein_term(p,o,e[i.I[1]],e[i.I[2]]), M))
end

function berstein_basis(prebasis,polytope::Polytope)
  C = _berstein_matrix(prebasis,polytope)
  linear_combination( C, prebasis )
end

berstein_basis(::MonomialBasis{0},::Polytope) = [ConstantField(1)]

function rationalize_bernstein_basis(basis,weights)
  ( basis .* weights ) ./ ( basis ⋅ weights )
end

