```@meta
CurrentModule = Gridap.ReferenceFEs
```

# Gridap.ReferenceFEs

```@docs
ReferenceFEs
``` 

## Polytopes

### Interface

```@docs
Polytope
get_faces(p::Polytope)
get_dimrange(p::Polytope)
Polytope{D}(p::Polytope,Dfaceid::Integer) where D
get_vertex_coordinates(p::Polytope)
(==)(a::Polytope{D},b::Polytope{D}) where D
get_edge_tangents(p::Polytope)
get_facet_normals(p::Polytope)
get_facet_orientations(p::Polytope)
get_vertex_permutations(p::Polytope)
get_vtkid(p::Polytope)
get_vtknodes(p::Polytope)
test_polytope(p::Polytope{D};optional::Bool) where D
num_dims(::Polytope)
num_faces(p::Polytope)
num_faces(p::Polytope,dim::Integer)
num_facets(p::Polytope)
num_edges(p::Polytope)
num_vertices(p::Polytope)
get_facedims(p::Polytope)
get_offsets(p::Polytope)
get_offset(p::Polytope,d::Integer)
get_faces(p::Polytope,dimfrom::Integer,dimto::Integer)
```
### Extrusion polytopes

```@docs
ExtrusionPolytope
ExtrusionPolytope(extrusion::Int...)
Polytope(extrusion::Int...)
HEX_AXIS
TET_AXIS
```

### Pre-defined polytope instances

```@docs
VERTEX
SEGMENT
TRI
QUAD
TET
HEX
WEDGE
PYRAMID
```
## Degrees of freedom

### Interface

```@docs
Dof
evaluate_dof!(cache,dof,field)
dof_cache(dof,field)
dof_return_type(dof,field)
test_dof(dof,field,v,comp::Function)
evaluate_dof(dof,field)
evaluate(dof::Dof,field)
```

### Working with arrays of DOFs

```@docs
evaluate_dof_array(dof::AbstractArray,field::AbstractArray)
evaluate(dof::AbstractArray{<:Dof},field::AbstractArray)
```
### Lagrangian dof bases

```@docs
LagrangianDofBasis
LagrangianDofBasis(::Type{T},nodes::Vector{<:Point}) where T
```

## Reference Finite Elements

### Interface

```@docs
ReferenceFE
num_dofs(reffe::ReferenceFE)
get_polytope(reffe::ReferenceFE)
get_prebasis(reffe::ReferenceFE)
get_dofs(reffe::ReferenceFE)
get_face_own_dofids(reffe::ReferenceFE)
ReferenceFE{N}(reffe::ReferenceFE,nfaceid::Integer) where N
get_own_dofs_permutations(reffe::ReferenceFE)
INVALID_PERM
get_shapefuns(reffe::ReferenceFE)
compute_shapefuns(dofs,prebasis)
num_dims(::ReferenceFE)
test_reference_fe(reffe::ReferenceFE{D}) where D
```

### Generic reference elements

```@docs
GenericRefFE
GenericRefFE(
  polytope::Polytope{D},
  prebasis::Field,
  dofs::Dof,
  facedofids::Vector{Vector{Int}};
  shapefuns::Field,
  ndofs::Int,
  dofperms::Vector{Vector{Int}},
  reffaces) where D
```

### Lagrangian reference elements

```@docs
LagrangianRefFE
LagrangianRefFE(
  polytope::Polytope{D},
  prebasis::MonomialBasis,
  dofs::LagrangianDofBasis,
  face_own_nodeids::Vector{Vector{Int}},
  own_nodes_permutations::Vector{Vector{Int}},
  reffaces...) where D
LagrangianRefFE(::Type{T},p::Polytope{D},orders) where {T,D}
compute_monomial_basis(::Type{T},p::Polytope,orders) where T
compute_own_nodes(p::Polytope,orders)
compute_face_orders(p::Polytope,face::Polytope,iface::Int,orders)
compute_nodes(p::Polytope,orders)
compute_own_nodes_permutations(p::Polytope, interior_nodes)
compute_lagrangian_reffaces(::Type{T},p::Polytope,orders) where T
num_nodes(reffe::LagrangianRefFE)
get_node_coordinates(reffe::LagrangianRefFE)
get_dof_to_node(reffe::LagrangianRefFE)
get_dof_to_comp(reffe::LagrangianRefFE)
get_node_and_comp_to_dof(reffe::LagrangianRefFE)
```
### Serendipity reference elements

```@docs
SerendipityRefFE
is_serendipity_compatible(p::Polytope)
```

