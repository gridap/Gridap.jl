
"""
    struct LagrangianRefFE{D} <: ReferenceFE{D}
      data::GenericRefFE{D}
      face_own_nodeids::Vector{Vector{Int}}
      own_nodes_permutations::Vector{Vector{Int}},
    end

Type representing a Lagrangian finite element. In addition to all the information
provided by a `ReferenceFE`, this type also provides "node-based" information in the following
fields:

- `face_own_nodeids::Vector{Vector{Int}}`: nodes owned by each face
- `own_nodes_permutations::Vector{Vector{Int}}`: permutations of the nodes when the vertices are permuted

For this type

-  `get_dofs(reffe)` returns a `LagrangianDofBasis`
-  `get_prebasis(reffe)` returns a `MonomialBasis`
-  `ReferenceFE{N}(reffe,faceid) where N` returns a `LagrangianRefFE{N}`

The following methods need to be overwritten in order to use the resulting reference FE in
a grid that is eventually written to vtk

- [`get_vtkid(p::Polytope,basis::MonomialBasis)`](@ref)
- [`get_vtknodes(p::Polytope,basis::MonomialBasis)`](@ref)

"""
struct LagrangianRefFE{D} <: ReferenceFE{D}
  data::GenericRefFE{D}
  face_own_nodeids::Vector{Vector{Int}}
  own_nodes_permutations::Vector{Vector{Int}}
  @doc """
      LagrangianRefFE(
        polytope::Polytope{D},
        prebasis::MonomialBasis,
        dofs::LagrangianDofBasis,
        face_own_nodeids::Vector{Vector{Int}},
        own_nodes_permutations::Vector{Vector{Int}},
        reffaces...) where D

  Low level (inner) constructor of `LagrangianRefFE`.
  """
  function LagrangianRefFE(
    polytope::Polytope{D},
    prebasis::MonomialBasis,
    dofs::LagrangianDofBasis,
    face_own_nodeids::Vector{Vector{Int}},
    own_nodes_permutations::Vector{Vector{Int}},
    reffaces...) where D

    face_own_dofids = _generate_nfacedofs(face_own_nodeids,dofs.node_and_comp_to_dof)
    own_dofs_permutations = _find_own_dof_permutaions(own_nodes_permutations,dofs.node_and_comp_to_dof,face_own_nodeids,face_own_dofids)

    data = GenericRefFE(
      polytope,prebasis,dofs,face_own_dofids;
      own_dofs_permutations = own_dofs_permutations,
      reffaces = reffaces)

    new{D}(data,face_own_nodeids,own_nodes_permutations)
  end
end

num_dofs(reffe::LagrangianRefFE) = reffe.data.ndofs

get_polytope(reffe::LagrangianRefFE) = reffe.data.polytope

get_prebasis(reffe::LagrangianRefFE) = reffe.data.prebasis

get_dofs(reffe::LagrangianRefFE) = reffe.data.dofs

get_face_own_dofids(reffe::LagrangianRefFE) = reffe.data.face_own_dofids

get_own_dofs_permutations(reffe::LagrangianRefFE) = reffe.data.own_dofs_permutations

get_shapefuns(reffe::LagrangianRefFE) = reffe.data.shapefuns

"""
    get_node_coordinates(reffe::LagrangianRefFE) -> Vector{Point{D,Float64}}

Returns the nodal coordinates of the underlying `LagrangianDofBasis`.
"""
get_node_coordinates(reffe::LagrangianRefFE) = reffe.data.dofs.nodes

"""
    get_dof_to_node(reffe::LagrangianRefFE) -> Vector{Int}

Returns the field `dof_to_node` of the underlying `LagrangianDofBasis`.
"""
get_dof_to_node(reffe::LagrangianRefFE) = reffe.data.dofs.dof_to_node

"""
    get_dof_to_comp(reffe::LagrangianRefFE) -> Vector{Int}

Returns the field `dof_to_comp` of the underlying `LagrangianDofBasis`.
"""
get_dof_to_comp(reffe::LagrangianRefFE) = reffe.data.dofs.dof_to_comp

"""
    get_node_and_comp_to_dof(reffe::LagrangianRefFE) -> Vector

Returns the field `node_and_comp_to_dof` of the underlying `LagrangianDofBasis`.
"""
get_node_and_comp_to_dof(reffe::LagrangianRefFE) = reffe.data.dofs.node_and_comp_to_dof

"""
    num_nodes(reffe::LagrangianRefFE) -> Int

Get the number of nodes in the Lagrangian reference FE
"""
num_nodes(reffe::LagrangianRefFE) = length(get_node_coordinates(reffe))

function ReferenceFE{N}(reffe::LagrangianRefFE,iface::Integer) where N
  ReferenceFE{N}(reffe.data,iface)
end

function ReferenceFE{D}(reffe::LagrangianRefFE{D},iface::Integer) where D
  @assert iface==1 "Only one D-face"
  reffe
end

# VTK related

"""
    get_vtkid(p::Polytope,basis::MonomialBasis) -> Int

Given a polytope `p` and a monomial basis, returns an integer with its vtk identifier.
Overloading of this function is needed only in order to visualize the underlying polytope
with Paraview.
"""
function get_vtkid(p::Polytope,basis::MonomialBasis)
  @abstractmethod
end

"""
    get_vtknodes(p::Polytope,basis::MonomialBasis) -> Vector{Int}

Given a polytope `p` and monomial basis, returns a vector of integers representing a permutation of the
polytope vertices required to relabel the vertices according the criterion adopted in
Paraview.
Overloading of this function is needed only in order to visualize the underlying polytope
with Paraview.
"""
function get_vtknodes(p::Polytope,basis::MonomialBasis)
  @abstractmethod
end

# Helpers for LagrangianRefFE

function _generate_nfacedofs(nfacenodes,node_and_comp_to_dof)
  faces = 1:length(nfacenodes)
  T = eltype(node_and_comp_to_dof)
  comps = 1:n_components(T)
  nfacedofs = [Int[] for i in faces]
  for face in faces
    nodes = nfacenodes[face]
    # Node major
    for comp in comps
      for node in nodes
        comp_to_dofs = node_and_comp_to_dof[node]
        dof = comp_to_dofs[comp]
        push!(nfacedofs[face],dof)
      end
    end
  end
  nfacedofs
end

function _find_own_dof_permutaions(node_perms,node_and_comp_to_dof,nfacenodeids,nfacedofsids)
  dof_perms = Vector{Int}[]
  T = eltype(node_and_comp_to_dof)
  ncomps = n_components(T)
  idof_to_dof = nfacedofsids[end]
  inode_to_node = nfacenodeids[end]
  for inode_to_pinode in node_perms
    ninodes = length(inode_to_pinode)
    nidofs = ncomps*ninodes
    idof_to_pidof = fill(INVALID_PERM,nidofs)
    for (inode,ipnode) in enumerate(inode_to_pinode)
      if ipnode == INVALID_PERM
        continue
      end
      node = inode_to_node[inode]
      pnode = inode_to_node[ipnode]
      comp_to_pdof = node_and_comp_to_dof[pnode]
      comp_to_dof = node_and_comp_to_dof[node]
      for comp in 1:ncomps
        dof = comp_to_dof[comp]
        pdof = comp_to_pdof[comp]
        idof = findfirst(i->i==dof,idof_to_dof)
        ipdof = findfirst(i->i==pdof,idof_to_dof)
        idof_to_pidof[idof] = ipdof
      end
    end
    push!(dof_perms,idof_to_pidof)
  end
  dof_perms
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
function LagrangianRefFE(::Type{T},p::Polytope{D},orders) where {T,D}
  prebasis = compute_monomial_basis(T,p,orders)
  nodes, face_own_nodeids = compute_nodes(p,orders)
  dofs = LagrangianDofBasis(T,nodes)
  interior_nodes = dofs.nodes[face_own_nodeids[end]]
  own_nodes_permutations = compute_own_nodes_permutations(p, interior_nodes)
  reffaces = compute_lagrangian_reffaces(T,p,orders)
  LagrangianRefFE(p,prebasis,dofs,face_own_nodeids,own_nodes_permutations,reffaces...)
end

function MonomialBasis(::Type{T},p::Polytope,orders) where T
  compute_monomial_basis(T,p,orders)
end

function LagrangianDofBasis(::Type{T},p::Polytope,orders) where T
  nodes, _ = compute_nodes(p,orders)
  LagrangianDofBasis(T,nodes)
end

# Constructors taking Int

function LagrangianRefFE(::Type{T},p::Polytope{D},order::Int) where {T,D}
  orders = tfill(order,Val{D}())
  LagrangianRefFE(T,p,orders)
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

    node_coords, face_own_nodeids = compute_nodes(p,orders)

Returns `node_coords`, the nodal coordinates of all the Lagrangian nodes associated with the order per direction
`orders`, and `face_own_nodeids`, being a vector of vectors indicating which nodes are owned by each of
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
  _compute_node_permutations(p, interior_nodes)
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
  if any( orders .== 0)
    _compute_constant_nodes(p,orders)
  elseif all(orders .== 1)
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
    face_shapefuns = change_basis(face_prebasis,change)
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
  lshapefuns = change_basis(lbasis,change)
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

function compute_monomial_basis(::Type{T},p::ExtrusionPolytope{D},orders) where {D,T}
  extrusion = Tuple(p.extrusion.array)
  terms = _monomial_terms(extrusion,orders)
  MonomialBasis{D}(T,orders,terms)
end

function compute_own_nodes(p::ExtrusionPolytope{D},orders) where D
  extrusion = Tuple(p.extrusion.array)
  if all(orders .== 0)
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
  nface = p.nfaces[iface+offset]
  face_orders = _extract_nonzeros(nface.extrusion,orders)
  face_orders
end

function compute_nodes(p::ExtrusionPolytope{D},orders) where D
  _nodes, facenodes = _compute_nodes(p,orders)
  if any( orders .== 0)
    return (_nodes, facenodes)
  end
  terms = _coords_to_terms(_nodes,orders)
  nodes = _terms_to_coords(terms,orders)
  (nodes, facenodes)
end

function get_vtkid(p::ExtrusionPolytope, basis::MonomialBasis)
  exponents = get_exponents(basis)
  vtkid, _ = _vtkinfo_extrusion_polytope(p,exponents)
  vtkid
end

function get_vtknodes(p::ExtrusionPolytope, basis::MonomialBasis)
  exponents = get_exponents(basis)
  _, vtknodes = _vtkinfo_extrusion_polytope(p,exponents)
  vtknodes
end

# Helpers for the ExtrusionPolytope-related implementation

function _monomial_terms(extrusion::NTuple{D,Int},orders) where D
  terms = CartesianIndex{D}[]
  if D == 0
    push!(terms,CartesianIndex(()))
    return terms
  end
  _check_orders(extrusion,orders)
  M = mutable(VectorValue{D,Int})
  term = zero(M)
  _orders = M(orders)
  k = 0
  _add_terms!(terms,term,extrusion,_orders,D,k)
  terms
end

function _interior_nodes(extrusion::NTuple{D,Int},orders) where D
  _check_orders(extrusion,orders)
  terms = CartesianIndex{D}[]
  M = mutable(VectorValue{D,Int})
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
  t = zero(mutable(P))
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
  x = zero(mutable(P))
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

function _vtkinfo_extrusion_polytope(p,exponents)

  n_nodes = length(exponents)

  if p == SEGMENT
    if n_nodes == 2
      vtkid = 3
      vtknodes = [1,2]
    else
      @notimplemented
    end

  elseif p == TRIANGLE
    if n_nodes == 3
      vtkid = 5
      vtknodes = [1,2,3]
    else
      @notimplemented
    end

  elseif p == QUAD
    if n_nodes == 4
      vtkid = 9
      vtknodes = [1,2,4,3]
    else
      @notimplemented
    end

  elseif p == TET
    if n_nodes == 4
      vtkid = 10
      vtknodes = [1,2,3,4]
    else
      @notimplemented
    end

  elseif p == HEX
    if n_nodes == 8
      vtkid = 12
      vtknodes = [1,2,4,3,5,6,8,7]
    else
      @notimplemented
    end

  else
    @notimplemented "vtkid not implemented for given ExtrusionPolytope"
  end

  (vtkid, vtknodes)
end

