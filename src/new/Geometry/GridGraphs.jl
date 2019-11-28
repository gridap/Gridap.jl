
"""
"""
abstract type GridGraph{D} end

## Discrete model
#
## DiscreteModel{Dc,Dp}
#
#num_faces(g,d)
#num_nodes(g)
#
#get_faces(g,dimfrom,dimto)
#
#get_vertex_node(g)
#get_node_vertex(g)
#
#get_face_nodes(g,d)#TODO harmonize with reffe
#get_faces_around_node(g,d)
#
#get_isboundary_face(g,d)
#get_isboundary_node(g)
#
#get_face_reffe_type(g,d) #TODO harmonize with grid
#get_face_polytope_type(g,d)
#
#get_reffes(ReferenceFE{d},g)
#get_polytopes(Polytope{d},g)
#
#get_vertex_coordinates(g)
#
#get_node_coordinates(g)
#
## Default implementations
#get_dface_to_faces(g,d)
#get_faces_around_dface(g,d)
#
#get_dimranges(g)
#get_offsets(g)
#get_offset(g,d)
#
#num_dims(g)
#num_dims(::Type)
#num_faces(g)
#num_cells(g)
#num_facets(g)
#num_edges(g)
#num_vertices(g)

"""
"""
function get_dimranges(g::GridGraph)
  @abstractmethod
end

"""
"""
function get_refcells(g::GridGraph)
  @abstractmethod
end

"""
"""
function get_cell_type(g::GridGraph)
  @abstractmethod
end

"""
"""
function get_cell_faces(g::GridGraph)
  @abstractmethod
end

"""
"""
function get_faces(g::GridGraph,dimfrom::Integer,dimto::Integer)
  @abstractmethod
end

"""
This one is optional and only needed to print the graph to vtk
"""
function get_vertex_coordinates(g::GridGraph)
  @abstractmethod
end

# Tests

function test_grid_graph(g::GridGraph;optional::Bool=false)
  D = num_dims(g)
  ranges = get_dimranges(g)
  @test length(ranges) == D+1
  @test isa(ranges,Vector{UnitRange{Int}})
  cell_to_faces = get_cell_faces(g)
  @test isa(cell_to_faces,AbstractArray{Vector{Int}})
  ncells = num_cells(g)
  @test length(cell_to_faces) == ncells
  nfaces = num_faces(g)
  cell_to_type = get_cell_type(g)
  @test isa(cell_to_type,AbstractVector{<:Integer})
  refcells = get_refcells(g)
  @test isa(refcells,AbstractVector{<:Polytope{D}})
  for dfrom in 0:D
    for dto in 0:D
       face_to_faces = get_faces(g,dfrom,dto)
       @test length(face_to_faces) == num_faces(g,dfrom)
       @test isa(face_to_faces,AbstractArray{Vector{Int}})
    end
  end
  if optional
    vertex_coods = get_vertex_coordinates(g)
    @test isa(vertex_coods,AbstractVector{<:Point{D}})
  end
end

# Some generic API

"""
"""
num_dims(g::GridGraph{D}) where D = D

num_dims(::Type{GridGraph{D}}) where D = D

"""
"""
num_cells(g::GridGraph) = length(get_cell_faces(g))

"""
"""
num_faces(g::GridGraph) = get_dimranges(g)[end][end]

"""
"""
num_faces(g::GridGraph,dim::Integer) = _num_faces(g,dim)

"""
"""
num_facets(g::GridGraph) = _num_facets(g)

"""
"""
num_edges(g::GridGraph) = _num_edges(g)

"""
"""
num_vertices(g::GridGraph) = _num_vertices(g)

"""
"""
get_facedims(g::GridGraph) = _get_facedims(g)

"""
"""
get_offsets(g::GridGraph) = _get_offsets(g)

"""
"""
get_offset(g::GridGraph,d::Integer) = _get_offset(g,d)

"""
"""
function is_boundary_face(g::GridGraph)
  @notimplemented
end

"""
"""
function is_boundary_face(g::GridGraph,facedim::Integer)
  @notimplemented
end

