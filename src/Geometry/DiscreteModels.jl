module DiscreteModels

using Test
using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery
using Gridap.GridGraphs: FullGridGraphFromData
using UnstructuredGrids: generate_dual_connections
using UnstructuredGrids: find_cell_to_faces

export DiscreteModel
export DiscreteModelFromFile
export test_discrete_model

import Gridap: Grid
import Gridap: GridGraph
import Gridap: FullGridGraph
import Gridap: FaceLabels
import Gridap: pointdim
import Gridap: celldim
import Gridap: Triangulation
import Gridap: labels_on_dim
import Gridap: tag_from_name
using JSON
import JSON: lower

# Interfaces

"""
D is number of components of the points in the model
"""
abstract type DiscreteModel{D} end

"""
extracts the Grid{D,Z} from the Model
"""
function Grid(::DiscreteModel{D},::Val{Z})::Grid{D,Z} where {D,Z}
  @abstractmethod
end

function GridGraph(::DiscreteModel)::GridGraph
  @abstractmethod
end

function FullGridGraph(::DiscreteModel)::FullGridGraph
  @abstractmethod
end

"""
Extracts the FaceLabels object providing information
about the geometrical labels and physical tags of all the
n-faces in the model for n=0,...,D
"""
function FaceLabels(::DiscreteModel{D})::FaceLabels{D} where D
  @abstractmethod
end

function labels_on_dim(model::DiscreteModel,dim::Integer)
  labels = FaceLabels(model)
  labels_on_dim(labels,dim)
end

Grid(m::DiscreteModel,dim::Integer) = Grid(m,Val(dim))

pointdim(::DiscreteModel{D}) where D = D

celldim(::DiscreteModel{D}) where D = D

function Triangulation(m::DiscreteModel,dim::Integer)
  grid = Grid(m,dim)
  Triangulation(grid)
end

function Triangulation(m::DiscreteModel{D}) where D
  Triangulation(m,D)
end

function tag_from_name(m::DiscreteModel,name::String)
  labels = FaceLabels(m)
  tag_from_name(labels,name)
end

Grid(model::DiscreteModel{D}) where D = Grid(model,D)

# Testers

function test_discrete_model(model::DiscreteModel{D},dim::Integer) where D
  @test dim == D
  @test pointdim(model) == D
  @test celldim(model) == D
  for d in 0:D
    grid = Grid(model,d)
    @test isa(grid,Grid)
    labels = FaceLabels(model)
    lab1 = labels_on_dim(model,d)
    lab2 = labels_on_dim(labels,d)
    @test lab1 == lab2
  end
  graph = GridGraph(model)
  fullgraph = FullGridGraph(model)
  @test isa(graph,GridGraph)
  @test isa(fullgraph,FullGridGraph)
  trian = Triangulation(model)
  @test isa(trian,Triangulation)
end

# Concrete implementations

struct DiscreteModelFromData{D} <: DiscreteModel{D}
  grids::Vector{Grid}
  graph::FullGridGraph
  facelabels::FaceLabels
end

function DiscreteModelFromData(
  grids::Vector{Grid},
  graph::FullGridGraph,
  facelabels::FaceLabels)

  D = ndims(graph) 
  @assert length(grids) == D + 1
  DiscreteModelFromData{D}(grids,graph,facelabels)
end

function Grid(model::DiscreteModelFromData{D},::Val{Z}) where {D,Z}
  model.grids[Z+1]
end

function GridGraph(model::DiscreteModelFromData)
  model.graph
end

function FullGridGraph(model::DiscreteModelFromData)
  model.graph
end

function FaceLabels(model::DiscreteModelFromData{D}) where D
  model.facelabels
end

# Serialization of discrete models

lower(model::DiscreteModel) = model_to_dict(model)

function model_to_dict(model::DiscreteModel{D}) where D

  dict = Dict{String,Any}()

  _add_dims!(dict,model)

  _add_nodes!(dict,model)

  for d in 1:D
    _add_grid!(dict,model,d)
  end

  facelabels = FaceLabels(model)

  _add_tags!(dict,facelabels)

  for d in 0:D
    _add_labels!(dict,facelabels,d)
  end

  dict

end

function _add_dims!(dict,model)
  D = celldim(model)
  Z = pointdim(model)
  dict["pointdim"] = Z
  dict["celldim"] = D
end

function _add_nodes!(dict,model)
  Z = pointdim(model)
  grid = Grid(model)
  node_to_coord = points(grid)
  nnodes = length(node_to_coord)
  node_to_coord_data = zeros(Float64,nnodes*Z)
  _fill_node_to_coord_data!(node_to_coord_data,node_to_coord)
  dict["nodes"] = node_to_coord_data
end

function _fill_node_to_coord_data!(node_to_coord_data,node_to_coord)
  i = 1
  for x in node_to_coord
    for xi in x
      node_to_coord_data[i] = xi
      i +=1
    end
  end
end

function _add_grid!(dict,model,d)

  grid = Grid(model,d)

  data, ptrs = compress(cells(grid))

  dict["face$(d)_data"] = data

  dict["face$(d)_ptrs"] = ptrs

  face_types, orders, extrusions =
    _setup_face_types(celltypes(grid),cellorders(grid))

  dict["face$(d)_types"] = face_types

  dict["orders$(d)"] = orders

  dict["extrusions$(d)"] = extrusions

end

function _setup_face_types(ct,co)
  @notimplemented
end

function _setup_face_types(ct::ConstantCellValue,co::ConstantCellValue)
  face_types = ones(Int,ct.length)
  (face_types, [co.value], [ct.value])
end

function _add_tags!(dict,facelabels)
  dict["tag_to_name"] = facelabels.tag_to_name
  dict["tag_to_labels"] = facelabels.tag_to_labels
end

function _add_labels!(dict,facelabels,d)
  face_to_label = labels_on_dim(facelabels,d)
  dict["labels$(d)"] = collect(face_to_label)
end

# Deserialization

function DiscreteModelFromFile(filename::AbstractString)

  base, extension = splitext(filename)
  s = Symbol(extension[2:end])
  DiscreteModelFromFile(filename,Val(s))

end

function DiscreteModelFromFile(filename::AbstractString,::Any)
  @notimplemented
end

function DiscreteModelFromFile(filename::AbstractString,::Val{:json})
    dict = JSON.parsefile(filename)
    model = dict_to_model(dict)
    model
end

function dict_to_model(dict::AbstractDict)

  facelabels = _setup_facelabels(dict)

  node_to_coord = _setup_node_to_coord(dict)

  grids = _setup_dim_to_grids(dict,node_to_coord)

  graph = _setup_graph(grids)

  DiscreteModelFromData(grids, graph, facelabels)

end

function _setup_facelabels(dict)

  D = dict["celldim"]

  dim_to_nface_to_label = Vector{Int}[]

  for d in 0:D
    a = dict["labels$(d)"]
    nnfaces = length(a)
    nface_to_label = zeros(Int,nnfaces)
    nface_to_label[:] = a
    push!(dim_to_nface_to_label,nface_to_label)
  end

  tag_to_name = convert(Vector{String},dict["tag_to_name"])

  tag_to_labels = convert(Vector{Vector{Int}},dict["tag_to_labels"])

  FaceLabels(
    dim_to_nface_to_label, tag_to_labels, tag_to_name)

end

function _setup_node_to_coord(dict)
  a = dict["nodes"]
  Z = dict["pointdim"]
  nnodes::Int = length(a)/Z
  T = Point{Z,Float64}
  node_to_coord = zeros(T,nnodes)
  m = zero(mutable(T))
  _fill_node_to_coord!(node_to_coord,a,m)
  node_to_coord
end

function _fill_node_to_coord!(node_to_coord,a,m)
  nnode = length(node_to_coord)
  i = 1
  for node in 1:nnode
    for comp in eachindex(m)
      m[comp] = a[i]
      i += 1
    end
    node_to_coord[node] = m
  end
end

function _setup_grid(dict,d,node_to_coord)
  data::Vector{Int} = dict["face$(d)_data"]
  ptrs::Vector{Int} = dict["face$(d)_ptrs"]
  orders = dict["orders$(d)"]
  extrusions = dict["extrusions$(d)"]
  @notimplementedif length(orders) != 1
  @notimplementedif length(extrusions) != 1
  extrusion::NTuple{d,Int} = tuple(extrusions[1]...)
  order::Int = orders[1]
  @notimplementedif order != 1 # We assume nodes == vertices in grid graph
  n = length(ptrs)-1
  ct = ConstantCellValue(extrusion,n)
  co = ConstantCellValue(order,n)
  UnstructuredGrid(node_to_coord,data,ptrs,ct,co)
end

function _setup_grid0(node_to_coord)
  n = length(node_to_coord)
  data = collect(1:n)
  ptrs = collect(1:(n+1))
  order = 1
  extrusion = ()
  ct = ConstantCellValue(extrusion,n)
  co = ConstantCellValue(order,n)
  UnstructuredGrid(node_to_coord,data,ptrs,ct,co)
end

function _setup_dim_to_grids(dict,node_to_coord)

  dim_to_grid = Grid[]

  grid0 = _setup_grid0(node_to_coord)
  push!(dim_to_grid,grid0)

  D = dict["celldim"]

  for d in 1:D

    grid = _setup_grid(dict,d,node_to_coord)

    push!(dim_to_grid,grid)

  end

  dim_to_grid
  
end

function _setup_graph(dim_to_grid)

  D = length(dim_to_grid)-1
  data = Matrix{IndexCellArray}(undef,D+1,D+1)

  for d in 0:D
    grid = dim_to_grid[d+1]
    nface_to_nodes = cells(grid)
    data[d+1,0+1] = nface_to_nodes
  end

  for d in 1:D
    nface_to_nodes = data[d+1,0+1]
    nfaces = length(nface_to_nodes)
    node_to_nfaces = generate_dual_connections(nface_to_nodes)
    data[0+1,d+1] = node_to_nfaces
    data[d+1,d+1] = _identity_cell_vector(nfaces)
  end

  for d in 1:D
    fgrid = dim_to_grid[d+1]
    for j in 1:(d-1)
      vertex_to_jfaces = data[0+1,j+1]
      face_to_jfaces = find_cell_to_faces(fgrid, vertex_to_jfaces, j)
      jface_to_faces = generate_dual_connections(face_to_jfaces)
      data[d+1,j+1] = face_to_jfaces
      data[j+1,d+1] = jface_to_faces
    end
  end

  FullGridGraphFromData(data)

end

function _identity_cell_vector(n)
  data = collect(1:n)
  ptrs = collect(1:(n+1))
  CellVectorFromDataAndPtrs(data,ptrs)
end

end # module
