
"""
    struct UnstructuredDiscreteModel{Dc,Dp,Tp,B} <: DiscreteModel{Dc,Dp}
      grid::UnstructuredGrid{Dc,Dp,Tp,B}
      grid_topology::UnstructuredGridTopology{Dc,Dp,Tp,B}
      face_labeling::FaceLabeling
    end
"""
struct UnstructuredDiscreteModel{Dc,Dp,Tp,B} <: DiscreteModel{Dc,Dp}
  grid::UnstructuredGrid{Dc,Dp,Tp,B}
  grid_topology::UnstructuredGridTopology{Dc,Dp,Tp,B}
  face_labeling::FaceLabeling
end

"""
    UnstructuredDiscreteModel(grid::Grid)
"""
function UnstructuredDiscreteModel(grid::Grid)
  _grid = UnstructuredGrid(grid)
  topo = UnstructuredGridTopology(_grid)
  nfaces = [num_faces(topo,d) for d in 0:num_cell_dims(topo)]
  labels = FaceLabeling(nfaces)
  UnstructuredDiscreteModel(_grid,topo,labels)
end

"""
"""
function UnstructuredDiscreteModel(model::DiscreteModel)
  grid = UnstructuredGrid(get_grid(model))
  topo = UnstructuredGridTopology(get_grid_topology(model))
  labels = get_face_labeling(model)
  UnstructuredDiscreteModel(grid,topo,labels)
end

function UnstructuredDiscreteModel(model::UnstructuredDiscreteModel)
  model
end

# Implementation of the interface

get_grid(model::UnstructuredDiscreteModel) = model.grid

get_grid_topology(model::UnstructuredDiscreteModel) = model.grid_topology

get_face_labeling(model::UnstructuredDiscreteModel) = model.face_labeling

# Io

function to_dict(model::UnstructuredDiscreteModel)
  dict = Dict{Symbol,Any}()
  grid = get_grid(model)
  labels = get_face_labeling(model)
  dict[:grid] = to_dict(grid)
  dict[:labeling] = to_dict(labels)
  dict[:vertex_node] = get_vertex_node(model)
  dict[:version] = "v0.7"
  D = num_cell_dims(model)
  dict[:D] = D
  topo = get_grid_topology(model)
  for d in 1:D
    k = Symbol("face_vertices_$d")
    dict[k] = to_dict(get_faces(topo,d,0))
  end
  dict
end

function from_dict(::Type{UnstructuredDiscreteModel},dict::Dict{Symbol,Any})

  if ! haskey(dict,:version)
    @unreachable "Cannot convert Dict to UnstructuredDiscreteModel for legacy format"
  end

  v = VersionNumber(dict[:version])
  if v < VersionNumber("v0.7")
    @unreachable "Cannot convert Dict to UnstructuredDiscreteModel for outdated format"
  end

  grid = from_dict(UnstructuredGrid,dict[:grid])
  labeling = from_dict(FaceLabeling,dict[:labeling])
  D::Int = dict[:D]
  vertex_to_node::Vector{Int} = dict[:vertex_node]
  orientation = OrientationStyle(grid)
  polytopes = map(get_polytope, get_reffes(grid))
  cell_type = get_cell_type(grid)
  vertex_coordinates = get_node_coordinates(grid)[vertex_to_node]

  nvertices = length(vertex_to_node)
  d_dface_to_vertices = Vector{Table{Int,Int32}}(undef,D+1)
  d_dface_to_vertices[0+1] = identity_table(Int,Int32,nvertices)
  for d in 1:D
    k = Symbol("face_vertices_$d")
    dface_to_vertices = from_dict(Table{Int,Int32},dict[k])
    d_dface_to_vertices[d+1] = dface_to_vertices
  end

  topo = UnstructuredGridTopology(
    vertex_coordinates,
    d_dface_to_vertices,
    cell_type,
    polytopes,
    orientation)

  UnstructuredDiscreteModel(grid,topo,labeling)
end

