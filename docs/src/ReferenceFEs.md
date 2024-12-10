```@meta
CurrentModule = Gridap.ReferenceFEs
```

# Gridap.ReferenceFEs

## Polytopes

```@docs
Polytope
Base.==
test_polytope
```

### Polytope Geometry

```@docs
get_edge_tangent
get_facet_normal
get_facet_orientations
get_vertex_coordinates
get_face_coordinates
get_bounding_box
```

### Polytope Topology

```@docs
num_dims
num_faces
num_facets
num_edges
num_vertices
get_faces
get_dimranges
get_dimrange
get_facedims
get_reffaces
get_face_dimranges
get_face_type
get_face_vertices
get_face_vertex_permutations
get_vertex_permutations
get_offsets
get_offset
is_n_cube
is_simplex
simplexify
```

### Extrusion Polytopes

```@docs
ExtrusionPolytope
get_extrusion
HEX_AXIS
TET_AXIS
VERTEX
SEGMENT
TRI
QUAD
TET
HEX
WEDGE
PYRAMID
next_corner!
```

### General Polytopes

```@docs
GeneralPolytope
Polygon
Polyhedron
get_graph
get_metadata
isopen
isactive
check_polytope_graph
simplexify_interior
simplexify_surface
```

## Quadratures

```@docs
Quadrature
get_coordinates
get_weights
get_name
num_points
num_point_dims
num_dims
test_quadratures
GenericQuadrature
```

### Available Quadratures

```@docs
TensorProduct
Duffy
Strang
XiaoGimbutas
```

## ReferenceFEs

```@docs
ReferenceFE
ReferenceFEName
GenericReferenceFE
num_dims
num_cell_dims
num_point_dims
num_faces
num_vertices
num_edges
num_facets
num_dofs
get_polytope
get_prebasis
get_shapefuns
get_dof_basis
compute_shapefuns
compute_dofs
Conformity
get_face_dofs
get_face_own_dofs
get_face_own_dofs_permutations
test_reference_fe
```

### Nodal ReferenceFEs

```@docs
LagrangianRefFE
is_first_order
is_P
is_Q
is_S
compute_monomial_basis
compute_own_nodes
compute_face_orders
compute_nodes
compute_own_nodes_permutations
compute_lagrangian_reffaces
VERTEX1
SEG2
QUAD4
TRI3
TET4
HEX8
```

### Moment-Based ReferenceFEs

```@docs
MomentBasedReferenceFE
RaviartThomasRefFE
NedelecRefFE
BDMRefFE
ArnoldWintherRefFE
MardalTaiWintherRefFE
```
