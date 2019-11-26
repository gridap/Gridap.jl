
abstract type GridGraph end

"""
"""
function num_dims(g::GridGraph)
  @abstractmethod
end

"""
"""
function get_dimranges(g::GridGraph)
  @abstractmethod
end

"""
"""
function get_cell_faces(g::GridGraph)
  @abstractmethod
end

"""
"""
function get_face_cells(g::GridGraph)
  @abstractmethod
end

"""
"""
function get_cell_faces(g::GridGraph,facedim::Integer)
  @abstractmethod
end

"""
"""
function get_face_cells(g::GridGraph,facedim::Integer)
  @abstractmethod
end

"""
"""
function get_vertex_node(g::GridGraph)
  @abstractmethod
end

# Tests

function test_grid_graph(g::GridGraph)
  D = num_dims(g)
  ranges = get_dimranges(g)
  @test length(ranges) == D+1
  @test isa(ranges,Vector{UnitRange{Int}})
  cell_to_faces = get_cell_faces(g)
  @test isa(cell_to_faces,AbstractArray{Vector{Int}})
  ncells = num_cells(g)
  @test length(cell_to_faces) == ncells
  face_to_cells = get_face_cells(g)
  @test isa(face_to_cells,AbstractArray{Vector{Int}})
  nfaces = num_faces(g)
  @test length(face_to_cells) == nfaces
  for d in 0:D
    cell_to_dfaces = get_cell_faces(g,d)
    @test isa(cell_to_dfaces,AbstractArray{Vector{Int}})
    dface_to_cells = get_face_cells(g,d)
    @test isa(dface_to_cells,AbstractArray{Vector{Int}})
  end
  vertex_to_node = get_vertex_node(g)
  @test isa(vertex_to_node,AbstractArray{<:Integer})
end

# Some generic API

"""
"""
num_cells(g::GridGraph) = length(get_cell_faces(g,0))

"""
"""
num_faces(g::GridGraph) = length(get_face_cells(g))

"""
"""
num_faces(g::GridGraph,dim::Integer) = GridGraph.ReferenceFEs._num_faces(g,dim)

"""
"""
num_facets(g::GridGraph) = Gridap.ReferenceFEs._num_facets(g)

"""
"""
num_edges(g::GridGraph) = Gridap.ReferenceFEs._num_edges(g)

"""
"""
num_vertices(g::GridGraph) = Gridap.ReferenceFEs._num_vertices(g)

"""
"""
get_facedims(g::GridGraph) = Gridap.ReferenceFEs._get_facedims(g)

"""
"""
get_offsets(g::GridGraph) = Gridap.ReferenceFEs._get_offsets(g)

"""
"""
get_offset(g::GridGraph,d::Integer) = Gridap.ReferenceFEs._get_offset(g,d)

