"""
    abstract type DiscreteModel{Dc,Dp} <: ConformingTriangulation{Dc,Dp}

Abstract type representing a conforming triangulation which
has information about the underlying geometical objects (i.e., the underlying
n-faces).

The `DiscreteModel` interfacy is defined by overloading the methods:

- [`get_faces(g::DiscreteModel,dimfrom::Integer,dimto::Integer)`](@ref)
- [`get_vertex_node(g::DiscreteModel)`](@ref)
- [`get_node_face_owner(g::DiscreteModel)`](@ref)
- [`get_face_nodes(g::DiscreteModel,d::Integer)`](@ref)
- [`get_isboundary_face(g::DiscreteModel,d::Integer)`](@ref)
- [`get_face_reffe_type(g::DiscreteModel,d::Integer)`](@ref)
- [`get_face_polytope_type(g::DiscreteModel,d::Integer)`](@ref)
- [`get_reffes(::Type{<:ReferenceFE{d}},g::DiscreteModel) where d`](@ref)
- [`get_polytopes(::Type{<:Polytope{d}},g::DiscreteModel) where d`](@ref)
- [`get_face_labeling(g::DiscreteModel)`](@ref)

and tested with this function:

- [`test_discrete_model`](@ref)

!!! warning
    Not sure if this abstract type is needed.
    We could include all the abstract methods here as optional abstract methods
    in the conforming triangulation interface. However, we have already
    introduced the concept of discrete model in the tutorials.
    This is a good argument to keep this type.

"""
abstract type DiscreteModel{Dc,Dp} <: ConformingTriangulation{Dc,Dp} end

"""
    get_faces(g::DiscreteModel,dimfrom::Integer,dimto::Integer)
"""
function get_faces(g::DiscreteModel,dimfrom::Integer,dimto::Integer)
  @abstractmethod
end

"""
    get_vertex_node(g::DiscreteModel)
"""
function get_vertex_node(g::DiscreteModel)
  @abstractmethod
end

"""
    get_node_face_owner(g::DiscreteModel)
"""
function get_node_face_owner(g::DiscreteModel)
  @abstractmethod
end

"""
    get_face_nodes(g::DiscreteModel,d::Integer)
"""
function get_face_nodes(g::DiscreteModel,d::Integer)
  @abstractmethod
end

"""
    get_face_own_nodes(g::DiscreteModel,d::Integer)
"""
function get_face_own_nodes(g::DiscreteModel,d::Integer)
  @abstractmethod
end

"""
    get_isboundary_face(g::DiscreteModel,d::Integer)
"""
function get_isboundary_face(g::DiscreteModel,d::Integer)
  @abstractmethod
end

"""
    get_face_reffe_type(g::DiscreteModel,d::Integer)

Index to the vector get_reffes(g,d)
"""
function get_face_reffe_type(g::DiscreteModel,d::Integer)
  @abstractmethod
end

"""
    get_face_polytope_type(g::DiscreteModel,d::Integer)
"""
function get_face_polytope_type(g::DiscreteModel,d::Integer)
  @abstractmethod
end

"""
    get_reffes(::Type{<:ReferenceFE{d}},g::DiscreteModel) where d
"""
function get_reffes(::Type{<:ReferenceFE{d}},g::DiscreteModel) where d
  @abstractmethod
end

"""
    get_polytopes(::Type{<:Polytope{d}},g::DiscreteModel) where d

Index to the vector get_polytopes(g,d)
"""
function get_polytopes(::Type{<:Polytope{d}},g::DiscreteModel) where d
  @abstractmethod
end

"""
    get_face_labeling(g::DiscreteModel)
"""
function get_face_labeling(g::DiscreteModel)
  @abstractmethod
end

# Implementation of ConformingTriangulation interface

function get_cell_nodes(g::DiscreteModel)
  D = num_cell_dims(g)
  get_face_nodes(g,D)
end

function get_reffes(g::DiscreteModel)
  D = num_cell_dims(g)
  get_reffes(NodalReferenceFE{D},g)
end

function get_cell_type(g::DiscreteModel)
  D = num_cell_dims(g)
  get_face_reffe_type(g,D)
end

# Default API

"""
    num_faces(g::DiscreteModel,d::Integer)
    num_faces(g::DiscreteModel)
"""
function num_faces(g::DiscreteModel,d::Integer)
  D = num_cell_dims(g)
  face_to_cells = get_faces(g,d,D)
  length(face_to_cells)
end

function num_faces(g::DiscreteModel)
  sum((num_faces(g,d) for d in 0:num_cell_dims(g)))
end

"""
    num_facets(g::DiscreteModel)
"""
num_facets(g::DiscreteModel) = _num_facets(g)

"""
    num_edges(g::DiscreteModel)
"""
num_edges(g::DiscreteModel) = _num_edges(g)

"""
    num_vertices(g::DiscreteModel)
"""
num_vertices(g::DiscreteModel) = _num_vertices(g)

"""
    get_vertex_coordinates(g::DiscreteModel)
"""
function get_vertex_coordinates(g::DiscreteModel)
  vertex_to_node = get_vertex_node(g)
  node_to_coords = get_node_coordinates(g)
  node_to_coords[vertex_to_node]
end

"""
    get_dimranges(g::DiscreteModel)
"""
function get_dimranges(g::DiscreteModel)
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
    get_offsets(g::DiscreteModel)
"""
get_offsets(g::DiscreteModel) = _get_offsets(g)

"""
    get_offset(g::DiscreteModel,d::Integer)
"""
get_offset(g::DiscreteModel,d::Integer) = _get_offset(g,d)

"""
    get_facedims(g::DiscreteModel)
"""
get_facedims(g::DiscreteModel) =  _get_facedims(Int8,g)

"""
    get_cell_faces(g::DiscreteModel)
"""
function get_cell_faces(g::DiscreteModel)
  D = num_cell_dims(g)
  d_to_cell_to_dface = Tuple((Table(get_faces(g,D,d)) for d in 0:D))
  offsets = Tuple(get_offsets(g))
  append_tables_locally(offsets,d_to_cell_to_dface)
end

"""
    get_isboundary_face(g::DiscreteModel)
"""
function get_isboundary_face(g::DiscreteModel)
  D = num_cell_dims(g)
  d_to_dface_to_isboundary = [get_isboundary_face(g,d) for d in 0:D]
  face_to_isboundary = vcat(d_to_dface_to_isboundary...)
end

"""
    get_isboundary_node(g::DiscreteModel)
"""
function get_isboundary_node(g::DiscreteModel)
  face_to_isboundary = get_isboundary_face(g)
  node_to_face_owner = get_node_face_owner(g)
  _get_isboundary_node(face_to_isboundary,node_to_face_owner)
end

function _get_isboundary_node(face_to_isboundary,node_to_face_owner)
  nnodes = length(node_to_face_owner)
  node_to_isboundary = fill(false,nnodes)
  for node in 1:nnodes
    face_owner = node_to_face_owner[node]
    isboundary = face_to_isboundary[face_owner]
    node_to_isboundary[node] = isboundary
  end
  node_to_isboundary
end

"""
    get_reffes_offsets(model::DiscreteModel)
"""
function get_reffes_offsets(model::DiscreteModel)
  D = num_cell_dims(model)
  d_to_offset = zeros(Int,D+1)
  for d in 1:D
    reffes = get_reffes(NodalReferenceFE{d-1},model)
    d_to_offset[d+1] = d_to_offset[d] + length(reffes)
  end
  d_to_offset
end

"""
    get_reffes_alldims(model::DiscreteModel)
"""
function get_reffes_alldims(model::DiscreteModel)
  D = num_cell_dims(model)
  vcat([get_reffes(NodalReferenceFE{d},model) for d in 0:D]...)
end

"""
    get_face_reffe_type(model::DiscreteModel)
"""
function get_face_reffe_type(model::DiscreteModel)
  D = num_cell_dims(model)
  data = [ get_face_reffe_type(model,d) for d in 0:D]
  offsets = get_reffes_offsets(model)
  for d in 1:D
    data[d+1] .+= offsets[d+1]
  end
  vcat(data...)
end

"""
    get_cell_perm_indices(model::DiscreteModel,d::Integer)
    get_cell_perm_indices(model::DiscreteModel)
"""
function get_cell_perm_indices(model::DiscreteModel)
  compute_cell_perm_indices(model)
end

function get_cell_perm_indices(model::DiscreteModel,d::Integer)
  compute_cell_perm_indices(model,d)
end

function compute_cell_perm_indices(model::DiscreteModel)
  D = num_cell_dims(model)
  data = [compute_cell_perm_indices(model,d) for d in 0:D]
  offsets = tfill(0,Val{D+1}()) 
  append_tables_locally(offsets,tuple(data...))
end

function compute_cell_perm_indices(model::DiscreteModel,d::Integer)

  D = num_cell_dims(model)
  face_to_fvertex_to_vertex = Table(get_faces(model,d,0))
  face_to_ftype = get_face_reffe_type(model,d)
  reffaces = get_reffes(NodalReferenceFE{d},model)
  facepolytopes = map(get_polytope,reffaces)
  ftype_to_pindex_to_cfvertex_to_fvertex = map(get_vertex_permutations,facepolytopes)
  cell_to_cvertex_to_vertex = Table(get_faces(model,D,0))
  cell_to_lface_to_face = Table(get_faces(model,D,d))
  cell_to_ctype = get_cell_type(model)
  reffes = get_reffes(model)
  polytopes = map(get_polytope,reffes)
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

function  _compute_cell_perm_indices!(
  cell_to_lface_to_pindex,
  cell_to_lface_to_face,
  cell_to_cvertex_to_vertex,
  cell_to_ctype,
  ctype_to_lface_to_cvertices,
  face_to_fvertex_to_vertex,
  face_to_ftype,
  ftype_to_pindex_to_cfvertex_to_fvertex)

  ncells = length(cell_to_lface_to_face)
  for cell in 1:ncells
    ctype = cell_to_ctype[cell]
    lface_to_cvertices = ctype_to_lface_to_cvertices[ctype]
    a = cell_to_lface_to_face.ptrs[cell]-1
    c = cell_to_cvertex_to_vertex.ptrs[cell]-1
    for (lface,cfvertex_to_cvertex) in enumerate(lface_to_cvertices)
      face = cell_to_lface_to_face.data[a+lface]
      ftype = face_to_ftype[face]
      b = face_to_fvertex_to_vertex.ptrs[face]-1
      pindex_to_cfvertex_to_fvertex = ftype_to_pindex_to_cfvertex_to_fvertex[ftype]
      pindexfound = false
      for (pindex, cfvertex_to_fvertex) in enumerate(pindex_to_cfvertex_to_fvertex)
        found = true
        for (cfvertex,fvertex) in enumerate(cfvertex_to_fvertex)
          vertex1 = face_to_fvertex_to_vertex.data[b+fvertex]
          cvertex = cfvertex_to_cvertex[cfvertex]
          vertex2 = cell_to_cvertex_to_vertex.data[c+cvertex]
          if vertex1 != vertex2
            found = false
            break
          end
        end
        if found
          cell_to_lface_to_pindex.data[a+lface] = pindex
          pindexfound = true
          break
        end
      end
      @assert pindexfound "Valid pindex not found"
    end
  end

end

function extract_face_reffes(
  ::Type{<:ReferenceFE{d}},
  model::DiscreteModel,
  reffes::Vector{<:NodalReferenceFE}) where d

  D = num_cell_dims(model)
  ctype_to_reffe = reffes
  ctype_to_lftype_to_refface = [
    get_reffes(NodalReferenceFE{d},reffe) for reffe in ctype_to_reffe]
  ctype_to_lface_to_lftype = [
    get_face_type(reffe,d) for reffe in ctype_to_reffe]
  t = _generate_ftype_to_refface(
    Val{d}(),ctype_to_lftype_to_refface,ctype_to_lface_to_lftype)
  ftype_to_refface, ctype_to_lface_to_ftype = t
  cell_to_faces = get_faces(model,D,d)
  cell_to_ctype = get_cell_type(model)
  nfaces = num_faces(model,d)

  face_to_ftype = generate_face_to_face_type(
    cell_to_faces, cell_to_ctype, ctype_to_lface_to_ftype, nfaces)

  (ftype_to_refface, face_to_ftype)
end

"""
    extract_face_reffes(
      ::Type{<:ReferenceFE{d}},
      model::DiscreteModel,
      reffes::Vector{<:NodalReferenceFE}) where d

    extract_face_reffes(
      model::DiscreteModel,reffes::Vector{<:NodalReferenceFE})
"""
function extract_face_reffes(
  model::DiscreteModel,reffes::Vector{<:NodalReferenceFE})
  D = num_cell_dims(model)
  facereffes = Vector{NodalReferenceFE}[]
  facetypes = Vector{Int8}[]
  for d in 0:D
    ftype_to_refface, face_to_ftype = extract_face_reffes(
      ReferenceFE{d},model,reffes)
    push!(facereffes,ftype_to_refface)
    push!(facetypes,face_to_ftype)
  end
  offsets = zeros(Int,D+1)
  for d in 1:D
    offsets[d+1] = offsets[d] + length(facereffes[d])
    facetypes[d+1] .+= offsets[d+1]
  end
  vcat(facereffes...), vcat(facetypes...)
end

"""
    test_discrete_model(model::DiscreteModel)
"""
function test_discrete_model(model::DiscreteModel{Dc,Dp}) where {Dc,Dp}
  test_conforming_triangulation(model)
  D = num_cell_dims(model)
  @test D == Dc
  @test num_point_dims(model) == Dp
  @test num_dims(model) == D
  vertex_to_node = get_vertex_node(model)
  @test isa(vertex_to_node,AbstractArray{<:Integer})
  @test length(vertex_to_node) == num_vertices(model)
  node_to_face_owner = get_node_face_owner(model)
  @test isa(node_to_face_owner,AbstractArray{<:Integer})
  @test length(node_to_face_owner) == num_nodes(model)
  labels = get_face_labeling(model)
  @test isa(labels,FaceLabeling)
  for n in 0:D
    nface_to_nodes = get_face_nodes(model,n)
    @test isa(nface_to_nodes,AbstractArray{<:Vector{<:Integer}})
    @test length(nface_to_nodes) == num_faces(model,n)
    nface_to_own_nodes = get_face_own_nodes(model,n)
    @test isa(nface_to_own_nodes,AbstractArray{<:Vector{<:Integer}})
    @test length(nface_to_own_nodes) == num_faces(model,n)
    nface_to_isboundary = get_isboundary_face(model,n)
    @test isa(nface_to_isboundary,AbstractArray{Bool})
    @test length(nface_to_isboundary) == num_faces(model,n)
    nface_to_ftype = get_face_reffe_type(model,n)
    @test isa(nface_to_ftype,AbstractArray{<:Integer})
    @test length(nface_to_ftype) == num_faces(model,n)
    nface_to_ftype = get_face_polytope_type(model,n)
    @test isa(nface_to_ftype,AbstractArray{<:Integer})
    @test length(nface_to_ftype) == num_faces(model,n)
    reffes = get_reffes(ReferenceFE{n},model)
    @test isa(reffes,Vector{<:ReferenceFE{n}})
    polytopes = get_polytopes(Polytope{n},model)
    @test isa(polytopes,Vector{<:Polytope{n}})
    for m in 0:D
      nface_to_mfaces = get_faces(model,n,m)
      @test isa(nface_to_mfaces,AbstractArray{<:Vector{<:Integer}})
      @test length(nface_to_mfaces) == num_faces(model,n)
    end
  end
end

# Low dim grids


"""
    ConformingTriangulation(::Type{<:ReferenceFE{d}},model::DiscreteModel) where d
"""
function ConformingTriangulation(::Type{<:ReferenceFE{d}},model::DiscreteModel) where d
  UnstructuredGrid(NodalReferenceFE{d},model)
end

"""
    UnstructuredGrid(::Type{<:ReferenceFE{d}},model::DiscreteModel) where d
"""
function UnstructuredGrid(::Type{<:ReferenceFE{d}},model::DiscreteModel) where d
  reffes = get_reffes(NodalReferenceFE{d},model)
  cell_type = get_face_reffe_type(model,d)
  cell_nodes = Table(get_face_nodes(model,d))
  node_coordinates = get_node_coordinates(model)
  UnstructuredGrid(
    node_coordinates,
    cell_nodes,
    reffes,
    cell_type)
end

"""
    replace_reffes(model::DiscreteModel,reffes::Vector{<:NodalReferenceFE})
"""
function replace_reffes(model::DiscreteModel,reffes::Vector{<:NodalReferenceFE})
  replace_reffes(UnstructuredDiscreteModel(model),reffes)
end


