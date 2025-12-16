
"""
    abstract type Polytope{D} <: GridapType

Abstract type representing a polytope (i.e., a polyhedron in arbitrary dimensions).

`D` is the polytope's dimension, defined as the smaller dimension of a flat
hypersurface containing the polytope, e.g. a triangle in a 3D space is a
`Polytope{2}` because it is contained in a 2D plane.
`D` differs from the embedding dimension, that is the number of coordinates of its vertices.

We adopt the [usual nomenclature](https://en.wikipedia.org/wiki/Polytope) for polytope-related objects.
All objects in a polytope (from vertices to the polytope itself) are called *n-faces* or simply *faces*.
The notation *n-faces* is used only when it is needed to refer to the object dimension n. Otherwise we simply
use *face*. In addition, we say

- vertex (pl. vertices): for 0-faces
- edge: for 1-faces
- facet: for (`D-1`)-faces

The following constants represent the main reference polytopes:
[`VERTEX`](@ref), [`SEGMENT`](@ref) (edge), [`TRI`](@ref) (triangle),
[`QUAD`](@ref) (quadrilateral), [`TET`](@ref) (tetrahedron), [`HEX`](@ref)
(hexahedron), [`WEDGE`](@ref) (triangular prism) and [`PYRAMID`](@ref).

The `Polytope` interface is defined by overloading the following functions

- [`get_faces(p::Polytope)`](@ref)
- [`get_dimranges(p::Polytope)`](@ref)
- [`Polytope{D}(p::Polytope,faceid::Integer) where D`](@ref)
- [`get_vertex_coordinates(p::Polytope)`](@ref)
- [`(==)(a::Polytope{D},b::Polytope{D}) where D`](@ref)

And optionally these ones:

- [`get_edge_tangent(p::Polytope)`](@ref)
- [`get_facet_normal(p::Polytope)`](@ref)
- [`get_facet_orientations(p::Polytope)`](@ref)
- [`get_vertex_permutations(p::Polytope)`](@ref)
- [`is_n_cube(p::Polytope)`](@ref)
- [`is_simplex(p::Polytope)`](@ref)
- [`simplexify(p::Polytope)`](@ref)

The interface can be tested with the function

- [`test_polytope`](@ref)

There are two implementations, [`ExtrusionPolytope`](@ref) for usual reference
polytopes ⊂ [0,1]ᴰ, and [`GeneralPolytope`](@ref).
"""
abstract type Polytope{D} <: GridapType end

# Mandatory

"""
    get_faces(p::Polytope) -> Vector{Vector{Int}}

Given a polytope `p` the function returns a vector of vectors
defining the *incidence* relation of the faces in the polytope.

Each face in the polytope receives a unique integer id. The id 1 is assigned
to the first 0-face. Consecutive increasing ids are assigned to the other
0-faces, then to 1-faces, and so on. The polytope itself receives the largest id
which coincides with `num_faces(p)`. For a face id `iface`, `get_faces(p)[iface]`
is a vector of face ids, corresponding to the faces that are *incident* with the face
labeled with `iface`. That is, faces that are either on its boundary or the face itself.
In this vector of incident face ids, faces are ordered by dimension, starting with 0-faces.
Within each dimension, the labels are ordered in a consistent way with the polyope object
for the face `iface` itself.

# Examples

```jldoctest
using Gridap.ReferenceFEs

faces = get_faces(SEGMENT)
println(faces)

# output
Array{Int,1}[[1], [2], [1, 2, 3]]
```

The constant [`SEGMENT`](@ref) is bound to a predefined instance of polytope
that represents a segment.
The face labels associated with a segment are `[1,2,3]`, being `1` and `2` for the vertices and
`3` for the segment itself. In this case, this function returns the vector of vectors
`[[1],[2],[1,2,3]]` meaning that vertex `1` is incident with vertex `1` (idem for vertex 2), and that
the segment (id `3`) is incident with the vertices `1` and `2` and the segment itself.

"""
function get_faces(p::Polytope)
  @abstractmethod
end

"""
    get_dimranges(p::Polytope) -> Vector{UnitRange{Int}}

Given a polytope `p` it returns a vector of ranges. The entry `d+1` in this vector
contains the range of face ids for the faces of dimension `d`.

# Examples

```jldoctest
using Gridap.ReferenceFEs

ranges = get_dimranges(SEGMENT)
println(ranges)

# output
UnitRange{Int}[1:2, 3:3]
```
Face ids for the vertices in the segment range from 1 to 2 (2 vertices),
the face ids for edges in the segment range from 3 to 3 (only one edge with id 3).

"""
function get_dimranges(p::Polytope)
  @abstractmethod
end

"""
    get_dimrange(p::Polytope, d::Integer)

Indices range of the `d`-dimensional faces of `p`. Equivalent to
[`get_dimranges(p)[d+1]`](@ref get_dimrange).
"""
function get_dimrange(p::Polytope, d::Integer)
  get_dimranges(p)[d+1]
end

"""
    Polytope{D}(p::Polytope, faceid::Integer)

Returns a `Polytope{D}` object representing the "reference" polytope of the
`D`-face of `p` with id `faceid`. The value `faceid` refers to the numeration
restricted to the dimension `D` (it starts with 1 for the first `D`-face).
"""
function Polytope{D}(p::Polytope, Dfaceid::Integer) where D
  @abstractmethod
end

"""
    get_vertex_coordinates(p::Polytope) -> Vector{Point{D,Float64}}

Return a vector of points containing the coordinates of the vertices of `p`.
"""
function get_vertex_coordinates(p::Polytope)
  @abstractmethod
end

"""
    (==)(a::Polytope{D},b::Polytope{D}) where D

Returns `true` if the polytopes `a` and `b` are equivalent. Otherwise, it
returns `false`.
Note that the operator `==` returns `false` by default for polytopes
of different dimensions. Thus, this function has to be overloaded only
for the case of polytopes `a` and `b` of same dimension.
"""
function (==)(a::Polytope{D}, b::Polytope{D}) where D
  @abstractmethod
end

function (==)(a::Polytope, b::Polytope)
  false
end

# Optional

"""
    get_edge_tangent(p::Polytope) -> Vector{VectorValue{D,Float64}}

Given a polytope `p`, returns a vector of `VectorValue` objects
representing the unit tangent vectors to the polytope edges, in the order of
[`get_dimrange(p, 1)`](@ref get_dimrange).
"""
function get_edge_tangent(p::Polytope)
  @abstractmethod
end

"""
    get_facet_normal(p::Polytope{D}) -> Vector{VectorValue{D,Float64}}

Return a vector of `VectorValue`s representing the unit *outward* normal vectors
to the facets of `p`.
"""
function get_facet_normal(p::Polytope)
  @abstractmethod
end

"""
    get_facet_orientations(p::Polytope) -> Vector{Int}

Return a vector of integers of length `num_facets(p)`.
The facets whose vertices are ordered consistently with the outwards normal
vector (w.r.t. right-hand rule convention) have value `1`, and the others `-1`.

See also [`get_facet_normal`](@ref) and [`get_faces`](@ref).
"""
function get_facet_orientations(p::Polytope)
  @abstractmethod
end

"""
    get_vertex_permutations(p::Polytope) -> Vector{Vector{Int}}

Returns a vector of vectors containing all admissible permutations of the
vertices of `p`. An admissible permutation is one such that, if the vertices of
the polytope are re-labeled according to this permutation, the resulting
polytope preserves the shape of the original one.

These permutation are used to iddentify every possile way a geometrical map may
permute the vertices of the boundary faces of a reference polytope (of dimension
≤ 3) into the physical one. Only the iddentity permutation `[1, 2, …, N]` is
returned for 3D polytopes with `N` vertices. Indeed, a (3≥D)-dimensional
physical polytope `𝓟` is only mapped by one physical map, but its faces are
mapped by the physical map of all adjascent elements.

# Examples

```jldoctest
using Gridap.ReferenceFEs

perms = get_vertex_permutations(SEGMENT)
println(perms)

# output
Array{Int,1}[[1, 2], [2, 1]]

perms = get_vertex_permutations(TET)
println(perms)

# output
Array{Int,1}[[1, 2], [2, 1]]
```
The first admissible permutation for a segment is `[1,2]`,i.e., the identity.
The second one is `[2,1]`, i.e., the first vertex is relabeled as `2` and the
second vertex is relabeled as `1`.

"""
function get_vertex_permutations(p::Polytope)
  @abstractmethod
end

"""
    is_simplex(p::Polytope) -> Bool
"""
function is_simplex(p::Polytope)
  @abstractmethod
end

"""
    is_n_cube(p::Polytope) -> Bool
"""
function is_n_cube(p::Polytope)
  @abstractmethod
end

"""
    simplexify(p::Polytope; kwargs...) -> ( Vector{Vector{Int}}, ref_simplex )

Returns a partition of `p` into simplices. The returned couple contains the
vector of each simplex connectivity array, and the reference simplex of the
partition simplices. The node coordinates are that of `p`.

`ExtrusionPolytope`s of dimension ≤3 support the `positive=true` kwarg.
If set to true, all resulting simplices will keep the orientation of the original polytope.
"""
function simplexify(p::Polytope; kwargs...)
  @abstractmethod
end

# Some generic API

"""
    num_point_dims(::Polytope)
    num_point_dims(::Type{<:Polytope})

Embedding dimension of a `Polytope` (type), i.e. number of coordinates of its vertices.
"""
num_point_dims(::Type{<:Polytope}) = @abstractmethod
num_point_dims(p::Polytope) = num_point_dims(typeof(p))

"""
    num_dims(::Type{<:Polytope{D}}) = D
    num_dims(p::Polytope{D}) = D

Returns `D`.
"""
num_dims(::Type{<:Polytope{D}}) where D = D
num_dims(p::Polytope) = num_dims(typeof(p))

"""
    num_cell_dims(::Type{<:Polytope{D}})
    num_cell_dims(p::Polytope{D})

Returns `D`.
"""
num_cell_dims(::Type{<:Polytope{D}}) where D = D
num_cell_dims(p::Polytope) = num_cell_dims(typeof(p))

"""
    num_faces(p::Polytope)

Returns the total number of faces of `p` (counting `p` itself).

# Examples

```jldoctest
num_faces(SEGMENT)

#output
3  # 2 vertices + 1 edge.
```
"""
function num_faces(p::Polytope)
  length(get_faces(p))
end

"""
    num_faces(p::Polytope, d::Integer)

Returns the number of faces of dimension `d` of `p`.
"""
function num_faces(p::Polytope, d::Integer)
  _num_faces(p, d)
end

function _num_faces(p, d)
  length(get_dimranges(p)[d+1])
end

"""
    num_facets(p::Polytope)

Returns the number of facets of `p`.
"""
function num_facets(p::Polytope)
  _num_facets(p)
end

function _num_facets(p)
  D = num_dims(p)
  if D > 0
    num_faces(p, D - 1)
  else
    0
  end
end

"""
    num_edges(p::Polytope)

Returns the number of edges of `p`.
"""
function num_edges(p::Polytope)
  _num_edges(p)
end

function _num_edges(p)
  D = num_dims(p)
  if D > 0
    num_faces(p, 1)
  else
    0
  end
end


"""
    num_vertices(p::Polytope)

Returns the number of vertices of `p`.
"""
function num_vertices(p::Polytope)
  _num_vertices(p)
end

function _num_vertices(p)
  num_faces(p, 0)
end

"""
    get_facedims(p::Polytope) -> Vector{Int}

Return a vector indicating the dimension of each face of `p`.

# Examples

```jldoctest
using Gridap.ReferenceFEs

dims = get_facedims(SEGMENT)
println(dims)

# output
[0, 0, 1]

```

The first two faces in the segment (the two vertices) have dimension 0 and the
third face (the segment itself) has dimension 1

"""
function get_facedims(p::Polytope)
  _get_facedims(Int, p)
end

function _get_facedims(::Type{T}, p) where T
  n = num_faces(p)
  facedims = zeros(T, n)
  dimrange = get_dimranges(p)
  for (i, r) in enumerate(dimrange)
    d = i - 1
    facedims[r] .= d
  end
  facedims
end

"""
    get_offsets(p::Polytope{D}) -> Vector{Int}

Return a vector whose d+1 entry is the integer offset between a face id in the
global numeration of `p`'s faces to the numeration of `p`'s faces of
dimension d, for d ∈ {0:`D`}.

# Examples

```jldoctest
using Gridap.ReferenceFEs

offsets = get_offsets(SEGMENT)
println(offsets)

# output
[0, 2]  # The first vertex is at 1+0, the first edge at 1+2, in get_faces(SEGMENT)
```
"""
function get_offsets(p::Polytope)
  _get_offsets(p)
end

function _get_offsets(p)
  D = num_dims(p)
  dimrange = get_dimranges(p)
  offsets = zeros(Int, D + 1)
  k = 0
  for i in 0:D
    d = i + 1
    offsets[d] = k
    r = dimrange[d]
    k += length(r)
  end
  offsets
end

"""
    get_offset(p::Polytope, d::Integer)

Return the integer offset between a face id in the global numeration of `p`'s
faces to the numeration of `p`'s faces of dimension `d`.

Equivalent to [`get_offsets(p)[d+1]`](@ref get_offsets).
"""
function get_offset(p::Polytope, d::Integer)
  _get_offset(p, d)
end

function _get_offset(p, d)
  get_offsets(p)[d+1]
end

"""
    get_faces(p::Polytope, dimfrom::Integer, dimto::Integer) -> Vector{Vector{Int}}

For `dimfrom >= dimto`, returns a vector that, for each face of `p` of dimension
`dimfrom`, stores a vector of the ids of faces of dimension `dimto` on its boundary.

For `dimfrom < dimto`, returns a vector that, for each face of `p` of dimension
`dimfrom`, stores a face ids vector containing the ids of the `dimto`-dimensional
faces adjascent to (touching) it.

The numerations used in this function are the ones restricted to each dimension.

# Examples

```jldoctest
using Gridap.ReferenceFEs

edge_to_vertices = get_faces(QUAD,1,0)
println(edge_to_vertices)

vertex_to_edges_around = get_faces(QUAD,0,1)
println(vertex_to_edges_around)

# output
Array{Int,1}[[1, 2], [3, 4], [1, 3], [2, 4]]
Array{Int,1}[[1, 3], [1, 4], [2, 3], [2, 4]]
```
"""
function get_faces(p::Polytope, dimfrom::Integer, dimto::Integer)
  if dimfrom >= dimto
    _get_faces_primal(p, dimfrom, dimto)
  else
    _get_faces_dual(p, dimfrom, dimto)
  end
end

function _get_faces_primal(p, dimfrom, dimto)
  dimrange = get_dimranges(p)
  r = dimrange[dimfrom+1]
  faces = get_faces(p)
  faces_dimfrom = faces[r]
  n = length(faces_dimfrom)
  faces_dimfrom_dimto = Vector{Vector{Int}}(undef, n)
  offset = get_offset(p, dimto)
  facefrom_dimranges = get_face_dimranges(p, dimfrom)
  for i in 1:n
    rto = facefrom_dimranges[i][dimto+1]
    faces_dimfrom_dimto[i] = faces_dimfrom[i][rto] .- offset
  end
  faces_dimfrom_dimto
end

function _get_faces_dual(p, dimfrom, dimto)
  tface_to_ffaces = get_faces(p, dimto, dimfrom)
  nffaces = num_faces(p, dimfrom)
  fface_to_tfaces = [Int[] for in in 1:nffaces]
  for (tface, ffaces) in enumerate(tface_to_ffaces)
    for fface in ffaces
      push!(fface_to_tfaces[fface], tface)
    end
  end
  fface_to_tfaces
end

"""
    get_face_dimranges(p::Polytope)
    get_face_dimranges(p::Polytope, d::Integer)

Return a vector containing, for each face of `p`, the result of
[`get_dimranges`](@ref) on the reference polytope of the face.

If `d` is given, the returned vector only contains the data for the
`d`-dimensional faces of `p`.

# Examples

2-faces of a pyramid

```jldoctest
get_face_dimranges(PYRAMID,2)

# output
5-element Vector{Vector{UnitRange{Int64}}}:
 [1:4, 5:8, 9:9]  # one QUAD 2-face
 [1:3, 4:6, 7:7]  # and four TRI 2-faces
 [1:3, 4:6, 7:7]
 [1:3, 4:6, 7:7]
 [1:3, 4:6, 7:7]
```
"""
function get_face_dimranges(p::Polytope, d::Integer)
  n = num_faces(p, d)
  rs = Vector{UnitRange{Int}}[]
  for i in 1:n
    f = Polytope{d}(p, i)
    r = get_dimranges(f)
    push!(rs, r)
  end
  rs
end

function get_face_dimranges(p::Polytope)
  rs = Vector{UnitRange{Int}}[]
  D = num_dims(p)
  for b in 0:D
    rs = vcat(rs, get_face_dimranges(p, d))
  end
  rs
end

"""
    get_face_vertices(p::Polytope) -> Vector{Vector{Int}}
    get_face_vertices(p::Polytope, d::Integer) -> Vector{Vector{Int}}

Return a vector containing, for each face of `p`, the vector of indices of the
vertices of that face.

If `d` is given, the returned vector only contains the data for the
`d`-dimensional faces of `p`.
"""
function get_face_vertices(p::Polytope, d::Integer)
  get_faces(p, d, 0)
end

function get_face_vertices(p::Polytope)
  face_vertices = Vector{Int}[]
  for d in 0:num_dims(p)
    dface_to_vertices = get_faces(p, d, 0)
    for vertices in dface_to_vertices
      push!(face_vertices, vertices)
    end
  end
  face_vertices
end

"""
    get_reffaces(p::Polytope) -> Vector{<:Polytope}
    get_reffaces(::Type{Polytope{d}}, p::Polytope) -> Vector{<:Polytope{d}}

Get a vector of the unique reference polytopes for the faces of `p`.

If `d` is given, the returned vector contains the polytopes for the
`d`-dimensional faces of `p` only.

# Examples

Get the unique polytopes constituting a triangle

```jldoctest
using Gridap.ReferenceFEs

reffaces = get_reffaces(TRI)

println(reffaces)

# output
ExtrusionPolytope{2}[VERTEX, SEGMENT, TRI]
```


Get the unique polytopes for the facets of a wedge.

```jldoctest
using Gridap.ReferenceFEs

reffaces = get_reffaces(Polytope{2},WEDGE)

println(reffaces)

# output
ExtrusionPolytope{2}[TRI, QUAD]
```
"""
function get_reffaces(::Type{Polytope{d}}, p::Polytope) where d
  ftype_to_refface, = _compute_reffaces_and_face_types(p, Val{d}())
  collect(ftype_to_refface)
end

function get_reffaces(p::Polytope)
  ftype_to_refface, = _compute_reffaces_and_face_types(p)
  collect(ftype_to_refface)
end

"""
    get_face_type(p::Polytope) -> Vector{Int}
    get_face_type(p::Polytope, d::Integer) -> Vector{Int}

Return a vector containing, for each face of `p`, the index in
`get_reffaces(p)` of the reference polytope of that face.

If `d` is given, the returned vector contains the indices for the
`d`-dimensional faces of `p` into `get_reffaces(Polytope{d},p)`.

# Examples

Get the unique polytopes for the facets of a wedge and identify of which
type each face is.

```jldoctest
using Gridap.ReferenceFEs

reffaces = get_reffaces(Polytope{2},WEDGE)

face_types = get_face_type(WEDGE,2)

println(reffaces)
println(face_types)

# output
Gridap.ReferenceFEs.ExtrusionPolytope{2}[TRI, QUAD]
[1, 1, 2, 2, 2]

```

The three first facets are of type `1`, i.e, `QUAD`, and the last ones of type `2`, i.e., `TRI`.
"""
function get_face_type(p::Polytope, d::Integer)
  _, iface_to_ftype = _compute_reffaces_and_face_types(p, Val{d}())
  iface_to_ftype
end

function get_face_type(p::Polytope)
  _, iface_to_ftype = _compute_reffaces_and_face_types(p)
  iface_to_ftype
end

function _compute_reffaces_and_face_types(p::Polytope, ::Val{d}) where d
  iface_to_refface = [Polytope{d}(p, iface) for iface in 1:num_faces(p, d)]
  _find_unique_with_indices(iface_to_refface)
end

function _compute_reffaces_and_face_types(p::Polytope)
  D = num_cell_dims(p)
  d_to_refdfaces = Vector{Polytope}[]
  d_to_dface_to_ftype = Vector{Int8}[]
  for d in 0:D
    reffaces, face_to_ftype = _compute_reffaces_and_face_types(p, Val(d))
    push!(d_to_refdfaces, reffaces)
    push!(d_to_dface_to_ftype, face_to_ftype)
  end
  d_to_offset = zeros(Int, D + 1)
  for d in 1:D
    d_to_offset[d+1] = d_to_offset[d] + length(d_to_refdfaces[d])
    d_to_dface_to_ftype[d+1] .+= d_to_offset[d+1]
  end
  (collect(vcat(d_to_refdfaces...)), vcat(d_to_dface_to_ftype...), d_to_offset)
end

function _find_unique_with_indices(a_to_b)
  T = eltype(a_to_b)
  u_to_b = T[]
  _find_unique!(u_to_b, a_to_b)
  a_to_u = zeros(Int, length(a_to_b))
  _find_indexin!(a_to_u, a_to_b, u_to_b)
  (u_to_b, a_to_u)
end

function _find_unique!(f::Vector, itr, pred::Function=(==))
  for i in itr
    found = false
    for fi in f
      if pred(i, fi)
        found = true
      end
    end
    if !found
      push!(f, i)
    end
  end
  f
end

function _find_indexin!(a_to_index, a_to_b, index_to_b, pred::Function=(==))
  for (a, b) in enumerate(a_to_b)
    for (index, _b) in enumerate(index_to_b)
      if pred(b, _b)
        a_to_index[a] = index
        break
      end
    end
  end
  a_to_index
end

"""
    get_bounding_box(p::Polytope{D})

Return a couple of `Point{D}`s defining a bounding box containing `p`. The box
is the `D`-cuboid with all edges parallel to the Cartesian axes, and whose
diametraly opposed vertices are the two returned vertices.
"""
function get_bounding_box(p::Polytope)
  vertex_to_coords = get_vertex_coordinates(p)
  get_bounding_box(vertex_to_coords)
end

get_bounding_box(points) = get_bounding_box(identity, points)

function get_bounding_box(f, points)
  P = typeof(f(first(points)))
  T = eltype(P)
  D = length(P)
  pmin = Point(tfill(T(Inf), Val{D}()))
  pmax = Point(tfill(T(-Inf), Val{D}()))
  for p in points
    fp = f(p)
    pmin = min.(pmin, fp)
    pmax = max.(pmax, fp)
  end
  (pmin, pmax)
end

"""
    get_measure(p::Polytope)
    get_measure(p::Polytope, vertex_coords)

Returns the measure of the polytope `p` with coordinates `vertex_coords`.
"""
function get_measure(p::Polytope)
  vertex_coords = get_vertex_coordinates(p)
  get_measure(p, vertex_coords)
end

function get_measure(p::Polytope, vertex_coords)
  conn, simplex = simplexify(p)
  m = eltype(eltype(coords))
  D = num_dims(p)
  VD = Val(D + 1)
  for ids in conn
    field = affine_map(ntuple(i -> vertex_coords[ids[i]], VD))
    m += meas(field.gradient)
  end
  m / factorial(D)
end

"""
    get_diameter(p::Polytope)
    get_diameter(p::Polytope, vertex_coords)
    get_diameter(vertex_coords)

Returns the diameter of the polytope `p` with coordinates `vertex_coords`, defined as
the maximum distance between any two vertices of the polytope.
"""
function get_diameter(p::Polytope)
  vertex_coords = get_vertex_coordinates(p)
  get_diameter(vertex_coords)
end

function get_diameter(p::Polytope, vertex_coords)
  get_diameter(vertex_coords)
end

function get_diameter(vertex_coords)
  h = zero(eltype(eltype(vertex_coords)))
  for p in vertex_coords
    for q in vertex_coords
      h = max(h, norm(p - q))
    end
  end
  h
end

"""
    get_face_vertex_permutations(p::Polytope)
    get_face_vertex_permutations(p::Polytope, d::Integer)

Return a vector containing, for each face of `p`, the result of
[`get_vertex_permutations`](@ref) on the reference polytope of the face.

If `d` is given, the returned vector only contains the data for the
`d`-dimensional faces of `p`.
"""
function get_face_vertex_permutations(p::Polytope, d::Integer)
  reffaces = [Polytope{d}(p, iface) for iface in 1:num_faces(p, d)]
  map(get_vertex_permutations, reffaces)
end

function get_face_vertex_permutations(p::Polytope)
  D = num_cell_dims(p)
  p = [get_face_vertex_permutations(p, d) for d in 0:D]
  vcat(p...)
end

"""
    get_face_coordinates(p::Polytope)
    get_face_coordinates(p::Polytope, d::Integer)


Return a vector containing, for each face of `p`, the coordinates of the
vertices of that face.

If `d` is given, the returned vector only contains the data for the
`d`-dimensional faces of `p`.
"""
function get_face_coordinates(p::Polytope, d::Integer)
  vert_to_coord = get_vertex_coordinates(p)
  face_to_vertices = get_faces(p, d, 0)
  collect(lazy_map(Broadcasting(Reindex(vert_to_coord)), face_to_vertices))
end

function get_face_coordinates(p::Polytope)
  D = num_cell_dims(p)
  v = Vector{Point{D,Float64}}[]
  for d in 0:D
    append!(v, get_face_coordinates(p, d))
  end
  v
end

# Aggregate own data into faces

function face_own_data_to_face_data(
  poly::Polytope{D}, face_own_data::AbstractVector{<:AbstractVector{T}}
) where {D,T}
  face_data = Vector{Vector{T}}(undef, num_faces(poly))
  for d in 0:D
    d_offset = get_offset(poly, d)
    for dface in 1:num_faces(poly, d)
      data = T[]
      for dd in 0:d
        dd_offset = get_offset(poly, dd)
        for ddface in get_faces(poly, d, dd)[dface]
          append!(data, face_own_data[ddface+dd_offset])
        end
      end
      face_data[dface+d_offset] = data
    end
  end
  return face_data
end

# Testers

"""
    test_polytope(p::Polytope{D}; optional::Bool=false) where D

Function that stresses out the functions in the `Polytope` interface.
It tests whether the function in the polytope interface are defined
for the given object, and whether they return objects of the expected type.
With `optional=false` (the default), only the mandatory functions are checked.
With `optional=true`, the optional functions are also tested.
"""
function test_polytope(p::Polytope{D}; optional::Bool=false) where D
  @test D == num_dims(p)
  faces = get_faces(p)
  @test isa(faces, Vector{Vector{Int}})
  @test num_faces(p) == length(faces)
  offsets = get_offsets(p)
  @test isa(offsets, Vector{Int})
  @test length(offsets) == D + 1
  dimrange = get_dimranges(p)
  @test isa(dimrange, Vector{UnitRange{Int}})
  @test length(dimrange) == D + 1
  @test p == p
  for d in 0:D
    for id in 1:num_faces(p, d)
      pd = Polytope{d}(p, id)
      @test isa(pd, Polytope{d})
    end
  end
  for dimfrom in 0:D
    for dimto in 0:D
      fs = get_faces(p, dimfrom, dimto)
      @test isa(fs, Vector{Vector{Int}})
    end
  end
  reffaces = get_reffaces(p)
  facetypes = get_face_type(p)
  @test isa(reffaces, Vector{<:Polytope})
  @test isa(facetypes, Vector{<:Integer})
  @test num_faces(p) == length(facetypes)
  x = get_vertex_coordinates(p)
  @test isa(x, Vector{Point{D,Float64}})
  @test length(x) == num_faces(p, 0)
  face_x = get_face_coordinates(p)
  @test isa(face_x, Vector{Vector{Point{D,Float64}}})
  @test length(face_x) == num_faces(p)
  if optional
    fn = get_facet_normal(p)
    @test isa(fn, Vector{VectorValue{D,Float64}})
    @test length(fn) == num_facets(p)
    or = get_facet_orientations(p)
    @test isa(or, Vector{Int})
    @test length(or) == num_facets(p)
    et = get_edge_tangent(p)
    @test isa(et, Vector{VectorValue{D,Float64}})
    @test length(et) == num_edges(p)
    @test isa(is_simplex(p), Bool)
    @test isa(is_n_cube(p), Bool)
  end
end
