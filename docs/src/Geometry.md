
```@meta
CurrentModule = Gridap.Geometry
```
# Gridap.Geometry

```@docs
Geometry
```

## Triangulations

### Interface

```@docs
Triangulation
get_cell_coordinates(trian::Triangulation)
get_reffes(trian::Triangulation)
get_cell_type(trian::Triangulation)
test_triangulation
num_cells(trian::Triangulation)
num_cell_dims(::Triangulation{Dc,Dp}) where {Dc,Dp}
num_point_dims(::Triangulation{Dc,Dp}) where {Dc,Dp}
num_dims(g::Triangulation{Dc}) where Dc
is_affine(trian::Triangulation)
has_straight_faces(trian::Triangulation)
get_cell_reffes(trian::Triangulation)
get_cell_shapefuns(trian::Triangulation)
get_cell_map(trian::Triangulation)
```
## ConformingTriangulations

### Interface

```@docs
ConformingTriangulation
OrientationStyle(::Type{<:ConformingTriangulation})
is_oriented(a::ConformingTriangulation)
ConformityStyle
RegularConformity
IrregularHConformity
IrregularPConformity
IrregularHPConformity
ConformityStyle(::Type{<:ConformingTriangulation})
get_node_coordinates(trian::ConformingTriangulation)
get_cell_nodes(trian::ConformingTriangulation)
test_conforming_triangulation
num_nodes(trian::ConformingTriangulation)
ConformingTriangulation(reffe::NodalReferenceFE)
ConformingTriangulation(T::Type{<:ReferenceFE{d}},p::Polytope) where d
ConformingTriangulation(T::Type{<:ReferenceFE{d}},trian::ConformingTriangulation) where d
replace_reffes(grid::ConformingTriangulation,reffes::Vector{<:NodalReferenceFE})
```

### UnstructuredGrids

```@docs
UnstructuredGrid
UnstructuredGrid(
  node_coordinates::Vector{Point{Dp,Tp}},
  cell_nodes::Table{Ti},
  reffes::Vector{<:NodalReferenceFE{Dc}},
  cell_types::Vector) where {Dc,Dp,Tp,Ti}
UnstructuredGrid(trian::ConformingTriangulation)
UnstructuredGrid(reffe::NodalReferenceFE)
UnstructuredGrid(::Type{<:ReferenceFE{d}},p::Polytope) where d
UnstructuredGrid(::Type{<:ReferenceFE{d}},trian::ConformingTriangulation) where d
UnstructuredGrid(x::AbstractArray{<:Point})
```

### CartesianGrids

```@docs
CartesianGrid
CartesianGrid(desc::CartesianDescriptor{D,T,F}) where {D,T,F}
CartesianGrid(domain,partition,map::Function)
get_cartesian_descriptor(a::CartesianGrid)
CartesianDescriptor
CartesianDescriptor(origin,sizes,partition,map::Function)
CartesianDescriptor(domain,partition,map::Function=identity)
```

## DiscreteModels

### Interface

```@docs
DiscreteModel
get_faces(g::DiscreteModel,dimfrom::Integer,dimto::Integer)
get_vertex_node(g::DiscreteModel)
get_node_face_owner(g::DiscreteModel)
get_face_nodes(g::DiscreteModel,d::Integer)
get_isboundary_face(g::DiscreteModel,d::Integer)
get_face_reffe_type(g::DiscreteModel,d::Integer)
get_face_polytope_type(g::DiscreteModel,d::Integer)
get_reffes(::Type{<:ReferenceFE{d}},g::DiscreteModel) where d
get_polytopes(::Type{<:Polytope{d}},g::DiscreteModel) where d
get_face_labeling(g::DiscreteModel)
num_faces(g::DiscreteModel,d::Integer)
num_facets(g::DiscreteModel)
num_edges(g::DiscreteModel)
num_vertices(g::DiscreteModel)
get_vertex_coordinates(g::DiscreteModel)
get_dimranges(g::DiscreteModel)
get_offsets(g::DiscreteModel)
get_offset(g::DiscreteModel,d::Integer)
get_facedims(g::DiscreteModel)
get_cell_faces(g::DiscreteModel)
get_isboundary_face(g::DiscreteModel)
get_isboundary_node(g::DiscreteModel)
get_reffes_offsets(model::DiscreteModel)
get_reffes_alldims(model::DiscreteModel)
get_face_reffe_type(model::DiscreteModel)
get_cell_perm_indices(model::DiscreteModel)
extract_face_reffes(model::DiscreteModel,reffes::Vector{<:NodalReferenceFE})
test_discrete_model
ConformingTriangulation(T::Type{<:ReferenceFE{d}},model::DiscreteModel) where d
UnstructuredGrid(T::Type{<:ReferenceFE{d}},model::DiscreteModel) where d
replace_reffes(model::DiscreteModel,reffes::Vector{<:NodalReferenceFE})
```

### FaceLabeling

```@docs
FaceLabeling
FaceLabeling(d_to_num_dfaces::Vector{Int})
num_dims(lab::FaceLabeling)
num_cell_dims(lab::FaceLabeling)
num_tags(lab::FaceLabeling)
num_entities(lab::FaceLabeling)
num_faces(lab::FaceLabeling,d::Integer)
num_faces(lab::FaceLabeling)
num_vertices(lab::FaceLabeling)
num_edges(lab::FaceLabeling)
num_facets(lab::FaceLabeling)
num_cells(lab::FaceLabeling)
get_face_entity(lab::FaceLabeling,d::Integer)
get_face_entity(lab::FaceLabeling)
get_tag_entities(lab::FaceLabeling,tag::Integer)
get_tag_entities(lab::FaceLabeling)
get_tag_name(lab::FaceLabeling,tag::Integer)
get_tag_name(lab::FaceLabeling)
get_tag_from_name(lab::FaceLabeling,name::String)
get_tag_from_name(lab::FaceLabeling)
add_tag!(lab::FaceLabeling,name::String,entities::Vector{<:Integer})
add_tag_from_tags!(lab::FaceLabeling, name::String, tags::Vector{Int})
```

### UnstructuredDiscreteModels

```@docs
UnstructuredDiscreteModel
UnstructuredDiscreteModel(trian::ConformingTriangulation)
```
### CartesianDiscreteModels

```@docs
CartesianDiscreteModel
CartesianDiscreteModel(desc::CartesianDescriptor{D,T,F}) where {D,T,F}
CartesianDiscreteModel(args...)
get_cartesian_descriptor(a::CartesianDiscreteModel)
get_polytope(::Type{<:Polytope{d}},model::CartesianDiscreteModel) where d
get_polytope(model::CartesianDiscreteModel)
```

