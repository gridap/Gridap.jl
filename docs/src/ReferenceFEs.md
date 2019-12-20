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
get_dimranges(p::Polytope)
get_dimrange(p::Polytope,d::Integer)
Polytope{D}(p::Polytope,Dfaceid::Integer) where D
get_vertex_coordinates(p::Polytope)
(==)(a::Polytope{D},b::Polytope{D}) where D
get_edge_tangents(p::Polytope)
get_facet_normals(p::Polytope)
get_facet_orientations(p::Polytope)
get_vertex_permutations(p::Polytope)
is_simplex(p::Polytope)
is_n_cube(p::Polytope)
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
get_face_vertices(p::Polytope,dim::Integer)
get_reffaces(::Type{<:Polytope{d}},p::Polytope) where d
get_face_type(p::Polytope,d::Integer)
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
get_dof_basis(reffe::ReferenceFE)
get_face_own_dofs(reffe::ReferenceFE)
get_face_own_dofs_permutations(reffe::ReferenceFE)
get_face_dofs(reffe::ReferenceFE)
INVALID_PERM
test_reference_fe(reffe::ReferenceFE{D}) where D
num_dims(reffe::ReferenceFE)
num_cell_dims(reffe::ReferenceFE)
num_point_dims(reffe::ReferenceFE)
num_faces(reffe::ReferenceFE)
num_vertices(reffe::ReferenceFE)
num_edges(reffe::ReferenceFE)
num_facets(reffe::ReferenceFE)
get_face_own_dofs(reffe::ReferenceFE,d::Integer)
get_face_own_dofs_permutations(reffe::ReferenceFE,d::Integer)
get_own_dofs_permutations(reffe::ReferenceFE)
get_face_dofs(reffe::ReferenceFE,d::Integer)
get_shapefuns(reffe::ReferenceFE)
compute_shapefuns(dofs,prebasis)
```

### Generic reference elements

```@docs
GenericRefFE
GenericRefFE(
  ndofs::Int,
  polytope::Polytope{D},
  prebasis::Field,
  dofs::Dof,
  face_own_dofs::Vector{Vector{Int}},
  face_own_dofs_permutations::Vector{Vector{Vector{Int}}},
  face_dofs::Vector{Vector{Int}},
  shapefuns::Field=compute_shapefuns(dofs,prebasis)) where D
```

## Node-based reference Finite Elements

### Interface

```@docs
NodalReferenceFE
get_node_coordinates(reffe::NodalReferenceFE)
get_node_and_comp_to_dof(reffe::NodalReferenceFE)
get_face_own_nodes(reffe::NodalReferenceFE)
get_face_own_nodes_permutations(reffe::NodalReferenceFE)
get_face_nodes(reffe::NodalReferenceFE)
test_nodal_reference_fe
num_nodes(reffe::NodalReferenceFE)
get_dof_to_node(reffe::NodalReferenceFE)
get_own_nodes_permutations(reffe::NodalReferenceFE)
get_vertex_node(reffe::NodalReferenceFE)
get_face_own_nodes(reffe::NodalReferenceFE,d::Integer)
get_face_own_nodes_permutations(reffe::NodalReferenceFE,d::Integer)
get_face_nodes(reffe::NodalReferenceFE,d::Integer)
```

### GenericNodalRefFE

```@docs
GenericNodalRefFE
```

### Lagrangian reference elements

```@docs
LagrangianRefFE
LagrangianRefFE(
  polytope::Polytope{D},
  prebasis::MonomialBasis,
  dofs::LagrangianDofBasis,
  face_own_nodes::Vector{Vector{Int}},
  own_nodes_permutations::Vector{Vector{Int}},
  reffaces) where D
LagrangianRefFE(::Type{T},p::Polytope{D},orders) where {T,D}
compute_monomial_basis(::Type{T},p::Polytope,orders) where T
compute_own_nodes(p::Polytope,orders)
compute_face_orders(p::Polytope,face::Polytope,iface::Int,orders)
compute_nodes(p::Polytope,orders)
compute_own_nodes_permutations(p::Polytope, interior_nodes)
compute_lagrangian_reffaces(::Type{T},p::Polytope,orders) where T
get_dof_to_comp(reffe::LagrangianRefFE)
ReferenceFE{N}(reffe::LagrangianRefFE,nfaceid::Integer) where N
(==)(a::LagrangianRefFE{D},b::LagrangianRefFE{D}) where D
get_reffaces(T::Type{<:ReferenceFE{d}},reffe::LagrangianRefFE) where d
get_face_type(reffe::LagrangianRefFE, d::Integer)
is_first_order(reffe::LagrangianRefFE)
is_affine(reffe::LagrangianRefFE)
```
### Serendipity reference elements

```@docs
SerendipityRefFE
```

### Pre-defined ReferenceFE instances

```@docs
SEG2
TRI3
QUAD4
TET4
HEX8
```
