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
    get_cell_to_faces(g::DiscreteModel)
"""
function get_cell_to_faces(g::DiscreteModel)
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
    if face_owner == UNSET
      @notimplemented
    end
    isboundary = face_to_isboundary[face_owner]
    node_to_isboundary[node] = isboundary
  end
  node_to_isboundary
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
  for n in 0:D
    nface_to_nodes = get_face_nodes(model,n)
    @test isa(nface_to_nodes,AbstractArray{<:Vector{<:Integer}})
    @test length(nface_to_nodes) == num_faces(model,n)
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

