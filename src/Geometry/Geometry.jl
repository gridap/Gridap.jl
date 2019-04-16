module Geometry

# Dependencies of this module

using Numa.Helpers
using Numa.FieldValues
using Numa.Polytopes
using Numa.RefFEs
using Numa.CellValues
using Numa.CellFunctions

# Functionality provided by this module

export Triangulation
export Grid
export GridGraph
export GridGraphFromData
export points
export cells
export triangulation
export celltypes
export cellorders
export celltovefs
export veftocells
export gridgraph
export geomap
export cellcoordinates
export cellbasis
export ncells

"""
Minimal interface for a mesh used for numerical integration
"""
abstract type Triangulation{Z,D} end

function cellcoordinates(::Triangulation{Z,D})::CellPoints{D} where {Z,D}
 @abstractmethod
end

# @fverdugo Return the encoded extrusion instead?
# If we restrict to polytopes that can be build form
# extrusion tuples, then we cannot accommodate polygonal elements, etc.
# This is a strong limitation IMHO that has to be worked out
"""
Returns the tuple uniquely identifying the Polytope of each cell
"""
function celltypes(::Triangulation{Z,D})::CellValue{NTuple{Z}} where {Z,D}
  @abstractmethod
end

cellorders(::Triangulation)::CellValue{Int} = @abstractmethod

function cellbasis(trian::Triangulation{Z,D}) where {Z,D}
  _cellbasis(trian,celltypes(trian),cellorders(trian))
end

function geomap(self::Triangulation)
  coords = cellcoordinates(self)
  basis = cellbasis(self)
  expand(basis,coords)
end

function ncells(self::Triangulation)
  coords = cellcoordinates(self)
  length(coords)
end

#@fverdugo make Z,D and D,Z consistent
"""
Abstract type representing a FE mesh a.k.a. grid
D is the dimension of the coordinates and Z is the dimension of the cells
"""
abstract type Grid{D,Z} end

function points(::Grid{D})::IndexCellValue{Point{D}} where D
  @abstractmethod
end

cells(::Grid)::IndexCellVector{Int} = @abstractmethod

triangulation(grid::Grid) = TriangulationFromGrid(grid)

"""
Abstract type that provides extended connectivity information associated with a grid.
This is the basic interface needed to distribute dof ids in
the construction of FE spaces.
"""
abstract type GridGraph end

celltovefs(::GridGraph)::IndexCellVector{Int} = @abstractmethod

veftocells(::GridGraph)::IndexCellVector{Int} = @abstractmethod

"""
Extracts the grid graph of the given grid
"""
gridgraph(::Grid)::GridGraph = @notimplemented

struct GridGraphFromData{C<:IndexCellVector{Int},V<:IndexCellVector{Int}} <: GridGraph
  celltovefs::C
  veftocells::V
end

celltovefs(self::GridGraphFromData) = self.celltovefs

veftocells(self::GridGraphFromData) = self.veftocells

# Submodules

include("Cartesian.jl")
include("Unstructured.jl")

# Helpers

_cellbasis( trian, ctypes, corders ) = @notimplemented

function _cellbasis(
  trian::Triangulation{Z,D},
  ctypes::ConstantCellValue{NTuple{Z,Int}},
  corders::ConstantCellValue{Int}) where {Z,D}
  ct = celldata(ctypes)
  co = celldata(corders)
  polytope = Polytope(Polytopes.PointInt{Z}(ct...))
  reffe = LagrangianRefFE{Z,ScalarValue}(polytope,fill(co,Z))
  basis = shfbasis(reffe)
  CellBasisFromSingleInterpolation(basis)
end

struct TriangulationFromGrid{D,Z,G<:Grid{D,Z}} <: Triangulation{Z,D}
  grid::G
end

function cellcoordinates(self::TriangulationFromGrid)
  CellVectorFromLocalToGlobal(cells(self.grid),points(self.grid))
end

celltypes(self::TriangulationFromGrid) = celltypes(self.grid)

cellorders(self::TriangulationFromGrid) = cellorders(self.grid)

end # module Geometry
