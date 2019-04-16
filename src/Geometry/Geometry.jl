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

# @santiagobadia : For sure, you can create your own polytope and fill all the
# info needed. It is a "static" polytope. We could express the current polytope
# as PolytopeByExtrusion and define the abstract polytope.
# In any case, I don't see what is the practical point for generating something
# like this. You will need to put a functional space on top of it and it does
# not work for a general polytope (in fact, it is even hard for pyramids, etc)
# I find it quite impractical.
# In any case, you could replace NTuple{Z} by e.g. P
# and in the implementation of a constructor, P = typeof(polytope), or an array
# of these values for a hybrid mesh... See also below (l.145).
"""
Returns the tuple uniquely identifying the Polytope of each cell
"""
function celltypes(::Triangulation{Z,D})::CellValue{NTuple{Z}} where {Z,D}
  @abstractmethod
end
# @santiagobadia : I think that the celltypes should be in the grid. The grid
# relies on it. E.g., the numbering being used in the cell vertices is a
# particular one, the one of the corresponding polytope numbering for its
# vertices. E.g., the Grid numbering being used in Fempar, Deal.ii, and Gid
# is different for quads.

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

# @santiagobadia : Do we want to call it vertices? I think it is much more
# meaningful than points... It is a point but it is more than that. It is the
# set of points that define the polytope as its convex hull...
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

# @santiagobadia : I would put this method in the interface of Grid...
"""
Extracts the grid graph of the given grid
"""
gridgraph(::Grid)::GridGraph = @notimplemented


struct GridGraphFromData{C<:IndexCellVector{Int},V<:IndexCellVector{Int}} <: GridGraph
  celltovefs::C
  veftocells::V
end

# @santiagobadia : I need to know whether the vef is a vertex, edge, or face, i.e., its dim.
# Do we want to provide a richer interface here? Or do we extract the polytope
# from the cell and use it.
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
  # @santiagobadia : For me the value of ctypes would be a Polytope, not
  # an NTuple... related to the comment above...
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

# @santiagobadia : In the future we will probably need a computecoordinates
# that takes a low order grid and increases order, using high-order mesh
# generation algorithms (untangling etc). Probably a too advanced topic yet...
function cellcoordinates(self::TriangulationFromGrid)
  CellVectorFromLocalToGlobal(cells(self.grid),points(self.grid))
end

celltypes(self::TriangulationFromGrid) = celltypes(self.grid)

cellorders(self::TriangulationFromGrid) = cellorders(self.grid)

end # module Geometry
