module Geometry

# Dependencies of this module

using Test
using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery: CellVectorFromLocalToGlobal
using Gridap.CellValuesGallery: CellValueFromArray

# Functionality provided by this module

export Grid
export test_grid
export points
export cells
export celltypes
export cellorders
export celltovefs
export veftocells
export geomap
export npoints
export FaceLabels
export labels_on_dim
export labels_on_tag
export ntags
export tag_from_name
export name_from_tag
export DiscreteModel
export celldim
export pointdim
export FullGridGraph
export connections
import Gridap: ncells
import Gridap: CellPoints
import Gridap: CellRefFEs
import Gridap: CellBasis
import Gridap: Triangulation

"""
Abstract type representing a FE mesh a.k.a. grid (i.e., there is a 1 to 1
correspondence between a grid and a conforming Lagrangian FE space with no
Dirichlet boundary conditions).
D is the dimension of the coordinates and Z is the dimension of the cells.
"""
abstract type Grid{D,Z} end

function points(::Grid{D})::IndexCellValue{<:Point{D}} where D
  @abstractmethod
end

cells(::Grid)::IndexCellArray{Int,1} = @abstractmethod

function celltypes(::Grid{D,Z})::CellValue{NTuple{Z,Int}} where {D,Z}
  @abstractmethod
end

cellorders(::Grid)::CellValue{Int} = @abstractmethod

celldim(::Grid{D,Z}) where {D,Z} = Z

pointdim(::Grid{D,Z}) where {D,Z} = D

ncells(g::Grid) = length(celltypes(g))

npoints(g::Grid) = length(points(g))

Triangulation(grid::Grid) = TriangulationFromGrid(grid)

function test_grid(grid::Grid{D,Z},np,nc) where {D,Z}
  @test isa(points(grid),IndexCellValue{<:Point{D}})
  @test isa(cells(grid),IndexCellVector{<:Integer})
  @test isa(celltypes(grid),CellValue{NTuple{Z,Int}})
  @test isa(cellorders(grid),IndexCellValue{Int})
  @test np == length(points(grid))
  @test nc == length(cells(grid))
  @test nc == length(celltypes(grid))
  @test nc == length(cellorders(grid))
end

"""
Classification of nfaces into geometrical and physical labels
"""
struct FaceLabels
  dim_to_nface_to_label::Vector{IndexCellValue}
  tag_to_labels::Vector{Vector{Int}}
  tag_to_name::Vector{String}
end

function FaceLabels(
  dim_to_nface_to_label::Vector{Vector{Int}},
  physlabel_to_labels::Vector{Vector{Int}},
  tag_to_name::Vector{String})

  cv = [ CellValueFromArray(v) for v in dim_to_nface_to_label ]
  FaceLabels(cv, physlabel_to_labels, tag_to_name)
end

"""
Returns an AbstractVector{Int} that represent the label for
each nface of dimension dim
"""
labels_on_dim(l::FaceLabels,dim::Integer) = l.dim_to_nface_to_label[dim+1]

"""
Returns a Vector{Int} with the labels associated with a given phystag
"""
labels_on_tag(l::FaceLabels,tag::Integer) = l.tag_to_labels[tag]

ntags(l::FaceLabels) = length(l.tag_to_labels)

function tag_from_name(l::FaceLabels,name::String)
  for tag in 1:ntags(l)
    if l.tag_to_name[tag] == name
      return tag
    end
  end
  0
end

name_from_tag(l::FaceLabels,tag::Integer) = l.tag_to_name[tag]

struct FullGridGraph
  data::Array{IndexCellArray,2}
end

function connections(g::FullGridGraph,from::Integer,to::Integer)
  # @assert from != to || (from == 0 && to == 0)
  g.data[from+1,to+1]
end

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

function FullGridGraph(::DiscreteModel)
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

function FaceLabels(model::DiscreteModel,dim::Integer)
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

# Helpers

_cellbasis( trian, ctypes, corders ) = @notimplemented

function _cellbasis(
  trian::Triangulation{Z,D},
  ctypes::ConstantCellValue{NTuple{Z,Int}},
  corders::ConstantCellValue{Int}) where {Z,D}
  # @santiagobadia : For me the value of ctypes would be a Polytope, not
  # an NTuple... related to the comment above...
  # @fverdugo yes. It will be refactorized accordingly
  ct = celldata(ctypes)
  co = celldata(corders)
  polytope = Polytope(ct...)
  reffe = LagrangianRefFE{Z,ScalarValue}(polytope,fill(co,Z))
  basis = shfbasis(reffe)
  ConstantCellMap(basis, length(ctypes)...)
end

struct TriangulationFromGrid{D,Z,G<:Grid{D,Z}} <: Triangulation{Z,D}
  grid::G
end

# @santiagobadia : In the future we will probably need a computecoordinates
# that takes a low order grid and increases order, using high-order mesh
# generation algorithms (untangling etc). Probably a too advanced topic yet...
# @fverdugo yes, sure. It will be just a factory function that returns objects
# that fit in the current interface. So the code is prepared for this extension.
function CellPoints(self::TriangulationFromGrid)
  CellVectorFromLocalToGlobal(cells(self.grid),points(self.grid))
end

function CellRefFEs(trian::TriangulationFromGrid)
  grid = trian.grid
  _build_cell_ref_fes(celltypes(grid),cellorders(grid))
end

function CellBasis(trian::TriangulationFromGrid)
  reffes = CellRefFEs(trian)
  _setup_cell_basis(reffes)
end

function _build_cell_ref_fes(codes,orders)
  @notimplemented
end

function _build_cell_ref_fes(
  codes::ConstantCellValue, orders::ConstantCellValue)
  code = codes.value
  D = length(code)
  order = orders.value
  polytope = Polytope(code)
  _orders = fill(order,D)
  reffe = LagrangianRefFE{D,Float64}(polytope,_orders)
  ConstantCellValue(reffe,length(codes))
end

function _setup_cell_basis(reffes)
  @notimplemented
end

function _setup_cell_basis(reffes::ConstantCellValue)
  reffe = reffes.value
  basis = shfbasis(reffe)
  ConstantCellMap(basis,reffes.length)
end

end # module Geometry
