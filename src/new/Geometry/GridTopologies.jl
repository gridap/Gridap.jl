
"""
    abstract type GridTopology{Dc,Dp}


Abstract type representing the topological information associated with a grid.

The `GridTopology` interface is defined by overloading the methods:

- [`get_faces(g::GridTopology,dimfrom::Integer,dimto::Integer)`](@ref)
- [`get_polytopes(g::GridTopology)`](@ref)
- [`get_cell_type(g::GridTopology)`](@ref)
- [`get_vertex_coordinates(g::GridTopology)`](@ref)

and tested with this function:

- [`test_grid_topology`](@ref)

"""
abstract type GridTopology{Dc,Dp} end

"""
    get_faces(g::GridTopology,dimfrom::Integer,dimto::Integer)
"""
function get_faces(g::GridTopology,dimfrom::Integer,dimto::Integer)
  @abstractmethod
end

"""
    get_polytopes(g::GridTopology)
"""
function get_polytopes(g::GridTopology)
  @abstractmethod
end

"""
    get_cell_type(g::GridTopology)
"""
function get_cell_type(g::GridTopology)
  @abstractmethod
end


"""
    get_vertex_coordinates(g::GridTopology)
"""
function get_vertex_coordinates(g::GridTopology)
  @abstractmethod
end

"""
    test_grid_topology(top::GridTopology)
"""
function test_grid_topology(top::GridTopology{Dc,Dp}) where {Dc,Dp}
  D = num_cell_dims(top)
  @test D == Dc
  @test num_point_dims(top) == Dp
  @test num_dims(top) == D
  cell_to_type = get_cell_type(top)
  @test isa(cell_to_type,AbstractVector{<:Integer})
  @test length(cell_to_type) == num_cells(top)
  polytopes = get_polytopes(top)
  @test isa(polytopes,Vector{<:Polytope{D}})
  get_isboundary_face(top)
  get_cell_faces(top)
  get_face_vertices(top)
  for n in 0:D
    compute_reffaces(Polytope{n},top)
    for m in 0:D
      nface_to_mfaces = get_faces(top,n,m)
      @test isa(nface_to_mfaces,AbstractArray{<:Vector{<:Integer}})
      @test length(nface_to_mfaces) == num_faces(top,n)
    end
  end
end

# Default API

"""
    num_cell_dims(::GridTopology) -> Int
    num_cell_dims(::Type{<:GridTopology}) -> Int
"""
num_cell_dims(::GridTopology{Dc,Dp}) where {Dc,Dp} = Dc
num_cell_dims(::Type{<:GridTopology{Dc,Dp}}) where {Dc,Dp} = Dc

"""
    num_point_dims(::GridTopology) -> Int
    num_point_dims(::Type{<:GridTopology}) -> Int
"""
num_point_dims(::GridTopology{Dc,Dp}) where {Dc,Dp} = Dp
num_point_dims(::Type{<:GridTopology{Dc,Dp}}) where {Dc,Dp} = Dp

"""
    num_dims(::GridTopology) -> Int
    num_dims(::Type{<:GridTopology}) -> Int

Equivalent to `num_cell_dims`.
"""
num_dims(g::GridTopology{Dc}) where Dc = Dc
num_dims(::Type{<:GridTopology{Dc}}) where Dc = Dc

"""
    num_faces(g::GridTopology,d::Integer)
    num_faces(g::GridTopology)
"""
function num_faces(g::GridTopology,d::Integer)
  D = num_cell_dims(g)
  face_to_cells = get_faces(g,d,D)
  length(face_to_cells)
end

function num_faces(g::GridTopology)
  sum((num_faces(g,d) for d in 0:num_cell_dims(g)))
end

"""
"""
num_cells(g::GridTopology) = num_faces(g,num_cell_dims(g))

"""
    num_facets(g::GridTopology)
"""
num_facets(g::GridTopology) = _num_facets(g)

"""
    num_edges(g::GridTopology)
"""
num_edges(g::GridTopology) = _num_edges(g)

"""
    num_vertices(g::GridTopology)
"""
num_vertices(g::GridTopology) = _num_vertices(g)

"""
    get_dimranges(g::GridTopology)
"""
function get_dimranges(g::GridTopology)
  ranges = UnitRange{Int}[]
  k = 1
  for d in 0:num_cell_dims(g)
    nf = num_faces(g,d)
    j = k+nf-1
    push!(ranges,k:j)
    k = j+1
  end
  ranges
end

"""
    get_dimrange(g::GridTopology,d::Integer)
"""
function get_dimrange(g::GridTopology,d::Integer)
  get_dimranges(g)[d+1]
end

"""
    get_offsets(g::GridTopology)
"""
get_offsets(g::GridTopology) = _get_offsets(g)

"""
    get_offset(g::GridTopology,d::Integer)
"""
get_offset(g::GridTopology,d::Integer) = _get_offset(g,d)

"""
    get_facedims(g::GridTopology)
"""
get_facedims(g::GridTopology) =  _get_facedims(Int8,g)

"""
    get_cell_faces(g::GridTopology)

Defaults to

    compute_cell_faces(g)
"""
function get_cell_faces(g::GridTopology)
  compute_cell_faces(g::GridTopology)
end

"""
    compute_cell_faces(g::GridTopology)
"""
function compute_cell_faces(g::GridTopology)
  D = num_cell_dims(g)
  d_to_cell_to_dface = Tuple((Table(get_faces(g,D,d)) for d in 0:D))
  offsets = Tuple(get_offsets(g))
  append_tables_locally(offsets,d_to_cell_to_dface)
end

"""
    get_face_vertices(g::GridTopology,d::Integer)
"""
function get_face_vertices(g::GridTopology,d::Integer)
  get_faces(g,d,0)
end

"""
    get_face_vertices(g::GridTopology)

Defaults to

    compute_face_vertices(g)
"""
function get_face_vertices(g::GridTopology)
  compute_face_vertices(g)
end

"""
    compute_face_vertices(g::GridTopology)
"""
function compute_face_vertices(g::GridTopology)
  D = num_cell_dims(g)
  d_to_dface_to_vertices = Tuple((Table(get_faces(g,d,0)) for d in 0:D))
  append_tables_globally(d_to_dface_to_vertices...)
end

"""
    get_cell_vertices(g::GridTopology)
"""
function get_cell_vertices(g::GridTopology)
  D = num_cell_dims(g)
  get_faces(g,D,0)
end

"""
    is_simplex(p::GridTopology) -> Bool
"""
function is_simplex(p::GridTopology)
  all(map(is_simplex, get_polytopes(p)))
end

"""
    is_n_cube(p::GridTopology) -> Bool
"""
function is_n_cube(p::GridTopology)
  all(map(is_n_cube, get_polytopes(p)))
end

"""
    get_reffaces(::Type{<:Polytope{d}}, g::GridTopology) where d

By default, it calls to `compute_reffaces`.
"""
function get_reffaces(::Type{<:Polytope{d}}, g::GridTopology) where d
  reffaces, _ = compute_reffaces(Polytope{d},g)
  reffaces
end

"""
    get_face_type(g::GridTopology,d::Integer)

By default, it calls to `compute_reffaces`.
"""
function get_face_type(g::GridTopology,d::Integer)
  _, face_to_ftype  = compute_reffaces(Polytope{d},g)
  face_to_ftype
end

"""
function compute_reffaces(::Type{<:Polytope{d}}, g::GridTopology) where d
"""
function compute_reffaces(::Type{<:Polytope{d}}, g::GridTopology) where d
  D = num_cell_dims(g)
  ctype_to_polytope = get_polytopes(g)
  ctype_to_lftype_to_refface = [ get_reffaces(Polytope{d},polytope) for polytope in ctype_to_polytope]
  ctype_to_lface_to_lftype = [ get_face_type(polytope,d) for polytope in ctype_to_polytope]
  t = _generate_ftype_to_refface(Val{d}(),ctype_to_lftype_to_refface,ctype_to_lface_to_lftype)
  ftype_to_refface, ctype_to_lface_to_ftype = t
  cell_to_faces = get_faces(g,D,d)
  cell_to_ctype = get_cell_type(g)
  nfaces = num_faces(g,d)
  face_to_ftype = generate_face_to_face_type(
    cell_to_faces, cell_to_ctype, ctype_to_lface_to_ftype, nfaces)

  (collect1d(ftype_to_refface), face_to_ftype)
end

function compute_reffaces(::Type{<:Polytope{D}}, g::GridTopology{D}) where D
  (get_polytopes(g), get_cell_type(g))
end

"""
    get_isboundary_face(g::GridTopology)
"""
function get_isboundary_face(g::GridTopology)
  compute_isboundary_face(g)
end

"""
    get_isboundary_face(g::GridTopology,d::Integer)
"""
function get_isboundary_face(g::GridTopology,d::Integer)
  compute_isboundary_face(g,d)
end

"""
    compute_isboundary_face(g::GridTopology)
"""
function compute_isboundary_face(g::GridTopology)
  D = num_cell_dims(g)
  d_to_dface_to_isboundary = [compute_isboundary_face(g,d) for d in 0:D]
  vcat(d_to_dface_to_isboundary...)
end

"""
    compute_isboundary_face(g::GridTopology,d::Integer)
"""
function compute_isboundary_face(g::GridTopology,d::Integer)
  D = num_cell_dims(g)
  if d == D-1
    _compute_isboundary_facet!(g)
  else
    _compute_isboundary_face!(g,d)
  end
end

function _compute_isboundary_facet!(g)
  D = num_cell_dims(g)
  d = D-1
  face_to_cells = get_faces(g,d,D)
  generate_facet_to_isboundary(face_to_cells)
end

function _compute_isboundary_face!(g,d)
  D = num_cell_dims(g)
  cell_to_faces = get_faces(g,D,d)
  cell_to_facets = get_faces(g,D,D-1)
  cell_to_ctype = get_cell_type(g)
  polytopes = get_polytopes(g)
  ctype_to_lface_to_lfacets = map( (p)->get_faces(p,d,D-1), polytopes)
  facet_to_isboundary = get_isboundary_face(g,D-1)
  nfaces = num_faces(g,d)

  generate_face_to_isboundary_from_cells(
    facet_to_isboundary,
    cell_to_faces,
    cell_to_facets,
    cell_to_ctype,
    ctype_to_lface_to_lfacets,
    nfaces)
end

"""
    get_cell_permutations(top::GridTopology)
"""
function get_cell_permutations(top::GridTopology)
  compute_cell_permutations(top)
end

"""
    get_cell_permutations(top::GridTopology,d::Integer)
"""
function get_cell_permutations(top::GridTopology,d::Integer)
  compute_cell_permutations(top,d)
end

"""
    compute_cell_permutations(top::GridTopology)
"""
function compute_cell_permutations(top::GridTopology)
  D = num_cell_dims(top)
  tables = (compute_cell_permutations(top,d) for d in 0:D)
  append_tables_locally(tables...)
end

"""
    compute_cell_permutations(top::GridTopology,d::Integer)
"""
function compute_cell_permutations(top::GridTopology,d::Integer)

  D = num_cell_dims(top)
  face_to_fvertex_to_vertex = Table(get_faces(top,d,0))
  face_to_ftype = get_face_type(top,d)
  reffaces = get_reffaces(Polytope{d},top)
  ftype_to_pindex_to_cfvertex_to_fvertex = map(get_vertex_permutations,reffaces)
  cell_to_cvertex_to_vertex = Table(get_faces(top,D,0))
  cell_to_lface_to_face = Table(get_faces(top,D,d))
  cell_to_ctype = get_cell_type(top)
  polytopes = get_polytopes(top)
  ctype_to_lface_to_cvertices = map( (p)->get_faces(p,d,0), polytopes )
  data = similar(cell_to_lface_to_face.data)
  ptrs = cell_to_lface_to_face.ptrs
  cell_to_lface_to_pindex = Table(data,ptrs)

  _compute_cell_perm_indices!(
    cell_to_lface_to_pindex,
    cell_to_lface_to_face,
    cell_to_cvertex_to_vertex,
    cell_to_ctype,
    ctype_to_lface_to_cvertices,
    face_to_fvertex_to_vertex,
    face_to_ftype,
    ftype_to_pindex_to_cfvertex_to_fvertex)

  cell_to_lface_to_pindex
end


