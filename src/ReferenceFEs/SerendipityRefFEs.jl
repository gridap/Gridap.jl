
"""
    struct Serendipity  <: ReferenceFEName
"""
struct Serendipity <: ReferenceFEName end

"""
    const serendipity = Serendipity()

Singleton of the [`Serendipity`](@ref) reference FE name.
"""
const serendipity = Serendipity()

"""
    SerendipityRefFE(::Type{T}, p::Polytope, order::Int)
    SerendipityRefFE(::Type{T}, p::Polytope, orders::Tuple)

Return a Lagrangian reference FE whose underlying approximation space is the
serendipity polynomial space 𝕊r of order `order`. Implemented on n-cubes with
homogneous order.

The type of the polytope `p` has to implement all the queries detailed in the
constructor [`LagrangianRefFE(::Type{T}, p::Polytope{D}, orders) where {T,D}`](@ref).

# Examples

```jldoctest
using Gridap.ReferenceFEs

order = 2
reffe = SerendipityRefFE(Float64,QUAD,order)

println( num_dofs(reffe) )

# output
8

```
"""
function SerendipityRefFE(::Type{T},p::Polytope,order::Int;
  poly_type=Monomial) where T

  @assert is_n_cube(p) "Polytope not compatible with serendipity elements"
  if order > 0
    sp = SerendipityPolytope(p)
  else
    sp = p
  end
  LagrangianRefFE(T,sp,order; poly_type)
end

function SerendipityRefFE(::Type{T},p::Polytope,orders::Tuple;
  poly_type=Monomial) where T

  order = first(orders)
  @assert all( orders .== order ) "Serendipity FEs must be isotropic, got orders $orders."
  SerendipityRefFE(T,p,order; poly_type)
end

function ReferenceFE(p::Polytope,::Serendipity,::Type{T},order;kwargs...) where T
  SerendipityRefFE(T,p,order;kwargs...)
end


# Helper private type
struct SerendipityPolytope{D,P} <: Polytope{D}
  hex::P
  SerendipityPolytope(p::Polytope{D}) where D = new{D,typeof(p)}(p)
end

# Implemented polytope interface

function get_faces(p::SerendipityPolytope)
  get_faces(p.hex)
end

function get_dimranges(p::SerendipityPolytope)
  get_dimranges(p.hex)
end

function Polytope{N}(p::SerendipityPolytope,Nfaceid::Integer) where N
  face_hex = Polytope{N}(p.hex, Nfaceid)
  SerendipityPolytope(face_hex)
end

function Polytope{D}(p::SerendipityPolytope{D},Dfaceid::Integer) where D
  @assert Dfaceid == 1 "Only one D-face"
  p
end

function get_vertex_coordinates(p::SerendipityPolytope)
  get_vertex_coordinates(p.hex)
end

function (==)(a::SerendipityPolytope{D},b::SerendipityPolytope{D}) where D
  a.hex == b.hex
end

function get_vertex_permutations(p::SerendipityPolytope)
  get_vertex_permutations(p.hex)
end

LagrangianRefFE(p::SerendipityPolytope) = LagrangianRefFE(p.hex)

is_simplex(p::SerendipityPolytope) = false

is_n_cube(p::SerendipityPolytope) = true

get_extrusion(p::SerendipityPolytope{D}) where D = Point(tfill(HEX_AXIS,Val{D}()))

# Implemented polytope interface for LagrangianRefFEs

function compute_monomial_basis(::Type{T},p::SerendipityPolytope{D},orders) where {T,D}
  MonomialBasis(Val(D),T,orders,_ser_filter)
end

function compute_poly_basis(::Type{T},p::SerendipityPolytope{D},orders,poly_type) where {T,D}
  FEEC_poly_basis
  r = iszero(D) ? 0 : first(orders)
  FEEC_poly_basis(Val(D),T,r,0,:S,poly_type) # SᵣΛ⁰(□ᴰ)
end

function compute_own_nodes(p::SerendipityPolytope{0},orders)
  compute_own_nodes(p.hex,orders)
end

function compute_own_nodes(p::SerendipityPolytope{1},orders)
  compute_own_nodes(p.hex,orders)
end

_own_s_filter(e,o) = ( sum( [ i for i in e ] ) <= o && all( [ i > 1 for i in e ] ) )

function _compute_own_s_nodes(orders)
  _terms = _define_terms(_q_filter,orders)
  _sort_by_nfaces!(_terms,orders)
  mask = _compute_filter_mask(_terms,_own_s_filter,orders)
  _terms = lazy_map(Reindex(_terms),mask)
  terms = map(t->CartesianIndex(Tuple(t)),_terms)
  _terms_to_coords(terms,orders)
end

function compute_own_nodes(p::SerendipityPolytope,orders)
  _compute_own_s_nodes(orders)
end

function compute_face_orders(
  p::SerendipityPolytope{D},face::SerendipityPolytope{N},iface::Int,orders) where {D,N}
  compute_face_orders(p.hex,face.hex,iface,orders)
end

