module Geometry

# Dependencies of this module

using Gridap
using Gridap.Helpers
using Gridap.CellValuesGallery: CellVectorFromLocalToGlobal
using Gridap.CellValuesGallery: CellValueFromArray

# Functionality provided by this module

export Grid
export GridGraph
export GridGraphFromData
export points
export cells
export celltypes
export cellorders
export celltovefs
export veftocells
export gridgraph
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

#@fverdugo make Z,D and D,Z consistent
"""
Abstract type representing a FE mesh a.k.a. grid
D is the dimension of the coordinates and Z is the dimension of the cells
"""
abstract type Grid{D,Z} end

# @santiagobadia : Do we want to call it vertices? I think it is much more
# meaningful than points... It is a point but it is more than that. It is the
# set of points that define the polytope as its convex hull...
# @fverdugo Just to clarify since (I don't know why) I removed the queries
# celltypes and cellorders from the Grid interface (I have added them again).
# Grid can have high order cells in the current design. Thus, vertices would not be a meaningful name...
# but of course we can change the name points (and also cells) if we find better names.
# At the beginning, I was thinking on an abstract type that represents only a linear Grid (i.e.,
# with polytope info but without cell order info), but
# I am not sure if it is needed since it is just a particular case of the current Grid...
# Moreover, for pure integer-based info, we have GridGraph
function points(::Grid{D})::IndexCellValue{Point{D}} where D
  @abstractmethod
end

cells(::Grid)::IndexCellArray{Int,1} = @abstractmethod

function celltypes(::Grid{D,Z})::CellValue{NTuple{Z}} where {D,Z}
  @abstractmethod
end

cellorders(::Grid)::CellValue{Int} = @abstractmethod

celldim(::Grid{D,Z}) where {D,Z} = Z

pointdim(::Grid{D,Z}) where {D,Z} = D

ncells(g::Grid) = length(celltypes(g))

npoints(g::Grid) = length(points(g))

Triangulation(grid::Grid) = TriangulationFromGrid(grid)

"""
Abstract type that provides extended connectivity information associated with a grid.
This is the basic interface needed to distribute dof ids in
the construction of FE spaces.
"""
abstract type GridGraph end

celltovefs(::GridGraph)::IndexCellArray{Int,1} = @abstractmethod

veftocells(::GridGraph)::IndexCellArray{Int,1} = @abstractmethod

# @santiagobadia : I would put this method in the interface of Grid...
# @fverdugo this would require define GridGraph before grid (which I find
# quite wird.) Anyway I find the current solution acceptable
# since this is julia and gridgraph is not a TBP of Grid...
"""
Extracts the grid graph of the given grid
"""
gridgraph(::Grid)::GridGraph = @notimplemented #@fverdugo Replace by GridGraph

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

Grid(m::DiscreteModel,dim::Integer) = Grid(m,Val(dim))

GridGraph(m::DiscreteModel,dim::Integer) = GridGraph(m,Val(dim))

pointdim(::DiscreteModel{D}) where D = D

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

#@fverdugo to be deleted together with (old) GridGraph
struct GridGraphFromData{C<:IndexCellArray{Int,1},V<:IndexCellArray{Int,1}} <: GridGraph
  celltovefs::C
  veftocells::V
end

# @santiagobadia : I need to know whether the vef is a vertex, edge, or face, i.e., its dim.
# Do we want to provide a richer interface here? Or do we extract the polytope
# from the cell and use it.
# @fverdugo yes sure. In fact I don't like to mix all vefs.
# I was thinking of an API like, e.g. for 3D,
# cell_to_faces = connections(graph,from=3,to=2)
# edge_to_cells = connections(graph,from=1,to=3).
# Once available, I would even delete celltovefs and veftocells since I don't like to mix things
# and I don't want to have duplicated data in some concrete implementations.
# (do you think it would be useful to keep them??)
# If we decide to keep them, I would propose an API like
# this one in order to be consistent.
# cell_to_vefs = connections(graph,from=3)
# vef_to_cells = connections(graph,to=3)
# moreover, we can also add
# face_to_mask = isonboundary(graph,dim=2)
# cell_to_mask = isonboundary(graph,dim=3)
celltovefs(self::GridGraphFromData) = self.celltovefs

veftocells(self::GridGraphFromData) = self.veftocells

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
  polytope = Polytope(Polytopes.PointInt{Z}(ct...))
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
