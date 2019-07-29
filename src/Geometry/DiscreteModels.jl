module DiscreteModels

using Test
using Gridap
using Gridap.Helpers

export DiscreteModel
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

end # module
