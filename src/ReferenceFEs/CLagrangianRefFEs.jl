struct GradConformity <: Conformity end
const H1Conformity = GradConformity


function Conformity(reffe::GenericLagrangianRefFE{GradConformity},sym::Symbol)
  h1 = (:H1,:C0,:Hgrad)
  if sym == :L2
    L2Conformity()
  elseif sym in h1
    H1Conformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a LagrangianRefFE with H1 conformity.

    Possible values of conformity for this reference fe are $((:L2, h1...)).
    """
  end
end

function Conformity(reffe::GenericLagrangianRefFE{L2Conformity},sym::Symbol)
  if sym == :L2
    L2Conformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a LagrangianRefFE with L2 conformity.

    Only conformity = :L2 allowed for this reference fe.
    """
  end
end

function get_face_own_nodes(reffe::GenericLagrangianRefFE{GradConformity},conf::GradConformity)
  p = get_polytope(reffe)
  orders = get_orders(reffe)
  nodes, face_own_nodes = compute_nodes(p,orders)
  face_own_nodes
end

function get_own_nodes_permutations(reffe::GenericLagrangianRefFE{GradConformity},conf::GradConformity)
  p = get_polytope(reffe)
  face_own_nodes = get_face_own_nodes(reffe)
  dofs = get_dof_basis(reffe)
  interior_nodes = dofs.nodes[face_own_nodes[end]]
  compute_own_nodes_permutations(p,interior_nodes)
end

function get_face_own_nodes_permutations(reffe::GenericLagrangianRefFE{GradConformity},conf::GradConformity)
  own_nodes_permutations = get_own_nodes_permutations(reffe)
  reffaces = reffe.reffe.metadata
  _reffaces = vcat(reffaces...)
  face_own_nodes_permutations = map(get_own_nodes_permutations,_reffaces)

  _compute_face_own_nodes_permutations(
    num_nodes(reffe),
    get_face_own_nodes(reffe),
    face_own_nodes_permutations,
    own_nodes_permutations)
end

function _compute_face_own_nodes_permutations(
  nnodes,
  face_own_nodes,
  face_own_nodes_permutations,
  own_nodes_permutations)

  if length(face_own_nodes_permutations) == 0
    # Vertex degenerated case
    return [own_nodes_permutations,]
  else
    if nnodes == length(face_own_nodes[end])
      # 0-order degenerated case
      _face_own_nodes_permutations = map( (x) -> fill(Int[],length(x)) , face_own_nodes_permutations )
    else
      # Standard case
      _face_own_nodes_permutations = copy(face_own_nodes_permutations)
    end
    push!(_face_own_nodes_permutations,own_nodes_permutations)
    return _face_own_nodes_permutations
  end
end

# API particular to LagrangianRefFE{GradConformity}

"""
    ReferenceFE{N}(reffe::GenericLagrangianRefFE{GradConformity},iface::Integer) where N
"""
function ReferenceFE{N}(reffe::GenericLagrangianRefFE{GradConformity},iface::Integer) where N
  reffaces = reffe.reffe.metadata
  reffaces[N+1][iface]
end

function ReferenceFE{D}(reffe::GenericLagrangianRefFE{GradConformity,D},iface::Integer) where D
  @assert iface==1 "Only one D-face"
  reffe
end

"""
    get_reffaces(
      ::Type{ReferenceFE{d}},
      reffe::GenericLagrangianRefFE{GradConformity}) where d -> Vector{GenericLagrangianRefFE{GradConformity,M,d}}
"""
function get_reffaces(::Type{ReferenceFE{d}},reffe::GenericLagrangianRefFE{GradConformity}) where d
  ftype_to_reffe, _ = _compute_reffes_and_face_types(reffe,Val{d}())
  ftype_to_reffe
  [reffe for reffe in ftype_to_reffe]
end

"""
    get_face_type(reffe::GenericLagrangianRefFE{GradConformity}, d::Integer) -> Vector{Int}
"""
function get_face_type(reffe::GenericLagrangianRefFE{GradConformity}, d::Integer)
  _, iface_to_ftype = _compute_reffes_and_face_types(reffe,Val{d}())
  iface_to_ftype
end

function _compute_reffes_and_face_types(reffe::GenericLagrangianRefFE{GradConformity},::Val{d}) where d
  p = get_polytope(reffe)
  iface_to_reffe = [ ReferenceFE{d}(reffe,iface) for iface in 1:num_faces(p,d) ]
  _find_unique_with_indices(iface_to_reffe)
end


"""
    is_first_order(reffe::GenericLagrangianRefFE{GradConformity}) -> Bool
"""
function is_first_order(reffe::GenericLagrangianRefFE{GradConformity})
  p = get_polytope(reffe)
  r = true
  r = r && num_vertices(p) == num_nodes(reffe)
  r = r && get_vertex_node(reffe) == collect(1:num_nodes(reffe))
  r
end

"""
    is_P(reffe::GenericLagrangianRefFE{GradConformity})
"""
function is_P(reffe::GenericLagrangianRefFE{GradConformity})
  is_simplex(get_polytope(reffe))
end

"""
   is_Q(reffe::GenericLagrangianRefFE{GradConformity})
"""
function is_Q(reffe::GenericLagrangianRefFE{GradConformity})
  monomials = get_prebasis(reffe)
  n = length(get_exponents(monomials))
  is_n_cube(get_polytope(reffe)) && (prod(get_orders(reffe).+1) == n)
end

"""
   is_S(reffe::GenericLagrangianRefFE{GradConformity})
"""
function is_S(reffe::GenericLagrangianRefFE{GradConformity})
  is_n_cube(get_polytope(reffe)) && ! is_Q(reffe)
end

function to_dict(reffe::GenericLagrangianRefFE{GradConformity})
  p = get_polytope(reffe)
  b = get_prebasis(reffe)
  dict = Dict{Symbol,Any}()
  dict[:orders] = collect(get_orders(reffe))
  dict[:extrusion] = Array(TensorValues.get_array(get_extrusion(p)))
  if is_S(reffe)
    dict[:space] = "serendipity"
  else
    dict[:space] = "default"
  end
  dict[:value] = string(return_type(b))
  dict
end

function from_dict(::Type{<:LagrangianRefFE},dict::Dict{Symbol,Any})
  orders = Tuple(dict[:orders])
  extrusion = Tuple(dict[:extrusion])
  if dict[:value] == "Float64"
    value = Float64
  else
    @notimplemented
  end
  space = dict[:space]
  p = Polytope(extrusion...)
  if space == "default"
    reffe = LagrangianRefFE(value,p,orders)
  elseif space == "serendipity"
    reffe = SerendipityRefFE(value,p,orders)
  else
    @unreachable "unknown space type"
  end
  reffe
end

# Construction of LagrangianRefFE from Polytopes

"""
    LagrangianRefFE(::Type{T},p::Polytope,orders) where T
    LagrangianRefFE(::Type{T},p::Polytope,order::Int) where T

Builds a `LagrangianRefFE` object on top of the given polytope. `T` is the type of
the value of the approximation space (e.g., `T=Float64` for scalar-valued problems,
`T=VectorValue{N,Float64}` for vector-valued problems with `N` components). The arguments `order` or `orders`
are for the polynomial order of the resulting space, which allows isotropic or anisotropic orders respectively
(provided that the cell topology allows the given anisotropic order). The argument `orders` should be an
indexable collection of `D` integers (e.g., a tuple or a vector), being `D` the number of space dimensions.

In order to be able to use this function, the type of the provided polytope `p` has to implement the
following additional methods. They have been implemented for `ExtrusionPolytope` in the library. They
need to be implemented for new polytope types in order to build Lagangian reference elements on top of them.

- [`compute_monomial_basis(::Type{T},p::Polytope,orders) where T`](@ref)
- [`compute_own_nodes(p::Polytope,orders)`](@ref)
- [`compute_face_orders(p::Polytope,face::Polytope,iface::Int,orders)`](@ref)

The following methods are also used in the construction of the `LagrangianRefFE` object. A default implementation
of them is available in terms of the three previous methods. However, the user can also implement them for
new polytope types increasing customization possibilities.

- [`compute_nodes(p::Polytope,orders)`](@ref)
- [`compute_own_nodes_permutations(p::Polytope, interior_nodes)`](@ref)
- [`compute_lagrangian_reffaces(::Type{T},p::Polytope,orders) where T`](@ref)
"""
function LagrangianRefFE(::Type{T},p::Polytope{D},orders;space::Symbol=_default_space(p)) where {T,D}
  if space == :P && is_n_cube(p)
    return _PDiscRefFE(T,p,orders)
  elseif space == :S && is_n_cube(p)
    SerendipityRefFE(T,p,orders)
  else
    if any(map(i->i==0,orders)) && !all(map(i->i==0,orders))
      cont = map(i -> i == 0 ? DISC : CONT,orders)
      return _cd_lagrangian_ref_fe(T,p,orders,cont)
    else
      return _lagrangian_ref_fe(T,p,orders)
    end
  end
end

function _default_space(p)
  if is_n_cube(p)
    :Q
  else
    :P
  end
end

function ReferenceFE(
  polytope::Polytope,
  ::Val{:Lagrangian},
  ::Type{T},
  orders::Union{Integer,Tuple{Vararg{Integer}}};
  space::Symbol=_default_space(polytope)) where T

  LagrangianRefFE(T,polytope,orders;space=space)
end


function _lagrangian_ref_fe(::Type{T},p::Polytope{D},orders) where {T,D}

  prebasis = compute_monomial_basis(T,p,orders)
  nodes, face_own_nodes = compute_nodes(p,orders)
  dofs = LagrangianDofBasis(T,nodes)
  reffaces = compute_lagrangian_reffaces(T,p,orders)

  nnodes = length(dofs.nodes)
  ndofs = length(dofs.dof_to_node)
  metadata = reffaces
  _reffaces = vcat(reffaces...)
  face_nodes = _generate_face_nodes(nnodes,face_own_nodes,p,_reffaces)
  face_own_dofs = _generate_face_own_dofs(face_own_nodes, dofs.node_and_comp_to_dof)
  face_dofs = _generate_face_dofs(ndofs,face_own_dofs,p,_reffaces)

  if all(map(i->i==0,orders) ) && D>0
    conf = L2Conformity()
  else
    conf = GradConformity()
  end

  reffe = GenericRefFE{typeof(conf)}(
    ndofs,
    p,
    prebasis,
    dofs,
    conf,
    metadata,
    face_dofs)

  GenericLagrangianRefFE(reffe,face_nodes)

end

function MonomialBasis(::Type{T},p::Polytope,orders) where T
  compute_monomial_basis(T,p,orders)
end

function LagrangianDofBasis(::Type{T},p::Polytope,orders) where T
  nodes, _ = compute_nodes(p,orders)
  LagrangianDofBasis(T,nodes)
end

# Helpers for LagrangianRefFE constructor

function _generate_face_nodes(nnodes,face_to_own_nodes,polytope,reffaces)

    face_to_num_fnodes = map(num_nodes,reffaces)
    push!(face_to_num_fnodes,nnodes)

    face_to_lface_to_own_fnodes = map(get_face_own_nodes,reffaces)
    push!(face_to_lface_to_own_fnodes,face_to_own_nodes)

    face_to_lface_to_face = get_faces(polytope)

  _generate_face_nodes_aux(
    nnodes,
    face_to_own_nodes,
    face_to_num_fnodes,
    face_to_lface_to_own_fnodes,
    face_to_lface_to_face)
end

function _generate_face_dofs(ndofs,face_to_own_dofs,polytope,reffaces)

    face_to_num_fdofs = map(num_dofs,reffaces)
    push!(face_to_num_fdofs,ndofs)

    face_to_lface_to_own_fdofs = map(get_face_own_dofs,reffaces)
    push!(face_to_lface_to_own_fdofs,face_to_own_dofs)

    face_to_lface_to_face = get_faces(polytope)

  _generate_face_nodes_aux(
    ndofs,
    face_to_own_dofs,
    face_to_num_fdofs,
    face_to_lface_to_own_fdofs,
    face_to_lface_to_face)
end

function _generate_face_nodes_aux(
  nnodes,
  face_to_own_nodes,
  face_to_num_fnodes,
  face_to_lface_to_own_fnodes,
  face_to_lface_to_face)

  if nnodes == length(face_to_own_nodes[end])
    face_fnode_to_node = fill(Int[],length(face_to_own_nodes))
    face_fnode_to_node[end] = collect(1:nnodes)
    return face_fnode_to_node
  end

  face_fnode_to_node = Vector{Int}[]
  for (face, nfnodes) in enumerate(face_to_num_fnodes)
    fnode_to_node = zeros(Int,nfnodes)
    lface_to_face = face_to_lface_to_face[face]
    lface_to_own_fnodes = face_to_lface_to_own_fnodes[face]
    for (lface, faceto) in enumerate(lface_to_face)
      own_nodes = face_to_own_nodes[faceto]
      own_fnodes = lface_to_own_fnodes[lface]
      fnode_to_node[own_fnodes] = own_nodes
    end
    push!(face_fnode_to_node,fnode_to_node)
  end

  face_fnode_to_node
end

# Constructors taking Int

function LagrangianRefFE(::Type{T},p::Polytope{D},order::Int;space::Symbol=_default_space(p)) where {T,D}
  orders = tfill(order,Val{D}())
  LagrangianRefFE(T,p,orders;space=space)
end

function MonomialBasis(::Type{T},p::Polytope{D},order::Int) where {D,T}
  orders = tfill(order,Val{D}())
  MonomialBasis(T,p,orders)
end

function LagrangianDofBasis(::Type{T},p::Polytope{D},order::Int) where {T,D}
  orders = tfill(order,Val{D}())
  LagrangianDofBasis(T,p,orders)
end

# Queries needed to be implemented for polytopes in order to use them
# for building LagrangianRefFEs in a seamless way

"""
    compute_monomial_basis(::Type{T},p::Polytope,orders) where T -> MonomialBasis

Returns the monomial basis of value type `T` and order per direction described by `orders`
on top of the polytope `p`.
"""
function compute_monomial_basis(::Type{T},p::Polytope,orders) where T
  @abstractmethod
end

"""
    compute_own_nodes(p::Polytope{D},orders) where D -> Vector{Point{D,Float64}}

Returns the coordinates of the nodes owned by the interior of the polytope
associated with a Lagrangian space with the order per direction described by `orders`.
"""
function compute_own_nodes(p::Polytope,orders)
  @abstractmethod
end

"""
    compute_face_orders(p::Polytope,face::Polytope,iface::Int,orders)

Returns a vector or a tuple with the order per direction at the face `face`
of the polytope `p` when restricting the order per direction `orders` to this face.
`iface` is the face id of `face` in the numeration restricted to the face dimension.
"""
function compute_face_orders(p::Polytope,face::Polytope,iface::Int,orders)
  @abstractmethod
end

"""
    compute_nodes(p::Polytope,orders)

When called

    node_coords, face_own_nodes = compute_nodes(p,orders)

Returns `node_coords`, the nodal coordinates of all the Lagrangian nodes associated with the order per direction
`orders`, and `face_own_nodes`, being a vector of vectors indicating which nodes are owned by each of
the faces of the polytope `p`.
"""
function compute_nodes(p::Polytope,orders)
  _compute_nodes(p,orders)
end

"""
    compute_own_nodes_permutations(
      p::Polytope, own_nodes_coordinates) -> Vector{Vector{Int}}

Returns a vector of vectors with the permutations of the nodes owned by the interior of the
polytope.
"""
function compute_own_nodes_permutations(p::Polytope, interior_nodes)
  perms = _compute_node_permutations(p, interior_nodes)
  perms
end

"""
    compute_lagrangian_reffaces(::Type{T},p::Polytope,orders) where T

Returns a tuple of length `D` being the number of space dimensions.
The entry `d+1` of this tuple contains a vector of `LagrangianRefFE`
one for each face of dimension `d` on the boundary of the polytope.
"""
function compute_lagrangian_reffaces(::Type{T},p::Polytope,orders) where T
  _compute_lagrangian_reffaces(T,p,orders)
end

# Default implementations

function _compute_nodes(p,orders)
  if any( map(i->i==0,orders))
    _compute_constant_nodes(p,orders)
  elseif all(map(i->i==1,orders))
    _compute_linear_nodes(p)
  else
    _compute_high_order_nodes(p,orders)
  end
end

function _compute_constant_nodes(p,orders)
  @assert all( orders .== 0) "If an order is 0 in some direction, it should be 0 also in the others"
  x = compute_own_nodes(p,orders)
  facenodes = [Int[] for i in 1:num_faces(p)]
  push!(facenodes[end],1)
  x, facenodes
end

function _compute_linear_nodes(p)
  x = get_vertex_coordinates(p)
  facenodes = [Int[] for i in 1:num_faces(p)]
  for i in 1:num_vertices(p)
    push!(facenodes[i],i)
  end
  x, facenodes
end

function _compute_high_order_nodes(p::Polytope{D},orders) where D
  nodes = Point{D,Float64}[]
  facenodes = [Int[] for i in 1:num_faces(p)]
  _compute_high_order_nodes_dim_0!(nodes,facenodes,p)
  for d in 1:(num_dims(p)-1)
    _compute_high_order_nodes_dim_d!(nodes,facenodes,p,orders,Val{d}())
  end
  _compute_high_order_nodes_dim_D!(nodes,facenodes,p,orders)
  (nodes, facenodes)
end

function _compute_high_order_nodes_dim_0!(nodes,facenodes,p)
  x = get_vertex_coordinates(p)
  k = 1
  for vertex in 1:num_vertices(p)
    push!(nodes,x[vertex])
    push!(facenodes[vertex],k)
    k += 1
  end
end

@noinline function _compute_high_order_nodes_dim_d!(nodes,facenodes,p,orders,::Val{d}) where d
  x = get_vertex_coordinates(p)
  offset = get_offset(p,d)
  k = length(nodes)+1
  for iface in 1:num_faces(p,d)
    face = Polytope{d}(p,iface)
    face_ref_x = get_vertex_coordinates(face)
    face_prebasis = MonomialBasis(Float64,face,1)
    change = inv(evaluate(face_prebasis,face_ref_x))
    face_shapefuns = linear_combination(change,face_prebasis)
    face_vertex_ids = get_faces(p,d,0)[iface]
    face_x = x[face_vertex_ids]
    face_orders = compute_face_orders(p,face,iface,orders)
    face_interior_nodes = compute_own_nodes(face,face_orders)
    face_high_x = evaluate(face_shapefuns,face_interior_nodes)*face_x
    for xi in 1:length(face_high_x)
      push!(nodes,face_high_x[xi])
      push!(facenodes[iface+offset],k)
      k += 1
    end
  end
end

function _compute_high_order_nodes_dim_D!(nodes,facenodes,p,orders)
  k = length(nodes)+1
  p_high_x = compute_own_nodes(p,orders)
  for xi in 1:length(p_high_x)
    push!(nodes,p_high_x[xi])
    push!(facenodes[end],k)
    k += 1
  end
end

_compute_node_permutations(::Polytope{0}, interior_nodes) = [[1]]

function _compute_node_permutations(p, interior_nodes)
  vertex_to_coord = get_vertex_coordinates(p)
  lbasis = MonomialBasis(Float64,p,1)
  change = inv(evaluate(lbasis,vertex_to_coord))
  lshapefuns = linear_combination(change,lbasis)
  perms = get_vertex_permutations(p)
  map = evaluate(lshapefuns,interior_nodes)
  pvertex_to_coord = similar(vertex_to_coord)
  node_perms = Vector{Int}[]
  tol = 1.0e-10
  for vertex_to_pvertex in perms
    node_to_pnode = fill(INVALID_PERM,length(interior_nodes))
    pvertex_to_coord[vertex_to_pvertex] = vertex_to_coord
    pinterior_nodes = map*pvertex_to_coord
    for node in 1:length(interior_nodes)
      x = interior_nodes[node]
      pnode = findfirst(i->norm(i-x)<tol,pinterior_nodes)
      if pnode != nothing
         node_to_pnode[node] = pnode
      end
    end
    push!(node_perms,node_to_pnode)
  end
  node_perms
end

_compute_lagrangian_reffaces(::Type{T},p::Polytope{0},orders) where T = ()

function _compute_lagrangian_reffaces(::Type{T},p::Polytope{D},orders) where {T,D}
  reffaces = [ LagrangianRefFE{d}[]  for d in 0:D ]
  p0 = Polytope{0}(p,1)
  reffe0 = LagrangianRefFE(T,p0,())
  for vertex in 1:num_vertices(p)
    push!(reffaces[0+1],reffe0)
  end
  offsets = get_offsets(p)
  for d in 1:(num_dims(p)-1)
    offset = offsets[d+1]
    for iface in 1:num_faces(p,d)
      face = Polytope{d}(p,iface)
      face_orders = compute_face_orders(p,face,iface,orders)
      refface = LagrangianRefFE(T,face,face_orders)
      push!(reffaces[d+1],refface)
    end
  end
  tuple(reffaces...)
end

# Particular implementation for ExtrusionPolytope

function LagrangianRefFE(p::ExtrusionPolytope)
  order = 1
  LagrangianRefFE(Float64,p,order)
end

function compute_monomial_basis(::Type{T},p::ExtrusionPolytope{D},orders) where {D,T}
  extrusion = Tuple(p.extrusion)
  terms = _monomial_terms(extrusion,orders)
  MonomialBasis{D}(T,orders,terms)
end

function compute_own_nodes(p::ExtrusionPolytope{D},orders) where D
  extrusion = Tuple(p.extrusion)
  if all(map(i->i==0,orders))
    _interior_nodes_order_0(p)
  else
    _interior_nodes(extrusion,orders)
  end
end

function _interior_nodes_order_0(p)
  x = get_vertex_coordinates(p)
  x0 = sum(x) / length(x)
  [x0,]
end

function compute_face_orders(p::ExtrusionPolytope,face::ExtrusionPolytope{D},iface::Int,orders) where D
  d = num_dims(face)
  offset = get_offset(p,d)
  nface = p.dface.nfaces[iface+offset]
  face_orders = _eliminate_zeros(Val{D}(),nface.extrusion,orders)
  Tuple(face_orders)
end

function _eliminate_zeros(::Val{d},a,o) where d
  b = zero(Mutable(Point{d,Int}))
  D = num_components(a)
  k = 1
  for i in 1:D
    m = a[i]
    if (m != 0)
      b[k] = o[i]
      k += 1
    end
  end
  Point(b)
end

function compute_nodes(p::ExtrusionPolytope{D},orders) where D
  _nodes, facenodes = _compute_nodes(p,orders)
  if any( map(i->i==0,orders))
    return (_nodes, facenodes)
  end
  terms = _coords_to_terms(_nodes,orders)
  nodes = _terms_to_coords(terms,orders)
  (nodes, facenodes)
end

# Helpers for the ExtrusionPolytope-related implementation

function _monomial_terms(extrusion::NTuple{D,Int},orders) where D
  terms = CartesianIndex{D}[]
  if D == 0
    push!(terms,CartesianIndex(()))
    return terms
  end
  _check_orders(extrusion,orders)
  M = Mutable(VectorValue{D,Int})
  term = zero(M)
  _orders = M(orders)
  k = 0
  _add_terms!(terms,term,extrusion,_orders,D,k)
  terms
end

function _interior_nodes(extrusion::NTuple{D,Int},orders) where D
  _check_orders(extrusion,orders)
  terms = CartesianIndex{D}[]
  M = Mutable(VectorValue{D,Int})
  term = zero(M)
  _orders = M(orders)
  k = 1
  _add_terms!(terms,term,extrusion,_orders,D,k)
  _terms_to_coords(terms,orders)
end

function _check_orders(extrusion,orders)
  D = length(extrusion)
  @assert length(orders) == D "container of orders not long enough"
  _orders = collect(orders)
  if extrusion[D] == HEX_AXIS
    _orders[D] = 0
  end
  for d in (D-1):-1:1
    if (extrusion[d] == HEX_AXIS || d == 1) && _orders[d+1] == 0
      _orders[d] = 0
    end
  end
  nz = _orders[_orders .!= 0]
  if length(nz) > 1
    @assert all(nz .== nz[1]) "The provided anisotropic order is not compatible with polytope topology"
  end
  nothing
end

function _add_terms!(terms,term,extrusion,orders,dim,k)
  _term = copy(term)
  _orders = copy(orders)
  indexbase = 1
  for i in k:(_orders[dim]-k)
    _term[dim] = i + indexbase
    if dim > 1
      if (extrusion[dim] == TET_AXIS) && i != 0
        _orders .-= 1
      end
      _add_terms!(terms,_term,extrusion,_orders,dim-1,k)
    else
      push!(terms,CartesianIndex(Tuple(_term)))
    end
  end
end

function _coords_to_terms(coords::Vector{<:Point{D}},orders) where D
  indexbase = 1
  terms = CartesianIndex{D}[]
  P = Point{D,Int}
  t = zero(Mutable(P))
  for x in coords
    for d in 1:D
      t[d] = round(x[d]*orders[d]) + indexbase
    end
    term = CartesianIndex(Tuple(t))
    push!(terms,term)
  end
  terms
end

function  _terms_to_coords(terms::Vector{CartesianIndex{D}},orders) where D
  P = Point{D,Float64}
  indexbase = 1
  nodes = P[]
  x = zero(Mutable(P))
  for t in terms
    for d in 1:D
      x[d] = (t[d] - indexbase) / orders[d]
    end
    node = P(x)
    push!(nodes,node)
  end
  nodes
end

function _extract_nonzeros(mask,values)
  b = Int[]
  for (m,n) in zip(mask,values)
    if (m != 0)
      push!(b, n)
    end
  end
  return Tuple(b)
end

# Precomputed instances

"""
    const VERTEX1 = LagrangianRefFE(Float64,VERTEX,1)
"""
const VERTEX1 = LagrangianRefFE(Float64,VERTEX,1)

"""
    const SEG2 = LagrangianRefFE(Float64,SEGMENT,1)
"""
const SEG2 = LagrangianRefFE(Float64,SEGMENT,1)

"""
    const QUAD4 = LagrangianRefFE(Float64,QUAD,1)
"""
const QUAD4 = LagrangianRefFE(Float64,QUAD,1)

"""
    const TRI3 = LagrangianRefFE(Float64,TRI,1)
"""
const TRI3 = LagrangianRefFE(Float64,TRI,1)

"""
    const TET4 = LagrangianRefFE(Float64,TET,1)
"""
const TET4 = LagrangianRefFE(Float64,TET,1)

"""
    const HEX8 = LagrangianRefFE(Float64,HEX,1)
"""
const HEX8 = LagrangianRefFE(Float64,HEX,1)
