module Geometry

# Dependencies of this module

using Numa.Helpers
using Numa.CellValues
using Numa.CellFunctions

# Functionality provided by this module

export IntegrationMesh
export Grid
export GridGraph
export GridGraphFromData
export points
export cells
export celltypes
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
abstract type IntegrationMesh{Z,D} end

function cellcoordinates(::IntegrationMesh{Z,D})::CellPoints{D} where {Z,D}
 @abstractmethod
end

function cellbasis(::IntegrationMesh{Z,D})::CellBasis{Z,Float64} where {Z,D}
  @abstractmethod
end

#@fverdugo TODO extend to high order
"""
Returns the tuple uniquely identifying the Polytope of each cell
"""
function celltypes(::IntegrationMesh{Z,D})::CellValue{NTuple{Z}} where {Z,D}
  @abstractmethod
end

function geomap(self::IntegrationMesh)
  coords = cellcoordinates(self)
  basis = cellbasis(self)
  expand(basis,coords)
end

function ncells(self::IntegrationMesh)
  coords = cellcoordinates(self)
  length(coords)
end


#@fverdugo TODO grid could also be high order
"""
Abstract type representing a FE mesh a.k.a. grid
D is the dimension of the coordinates and Z is the dimension of the cells
"""
abstract type Grid{D,Z} end

function points(::Grid{D})::IndexCellValue{Point{D}} where D
  @abstractmethod
end

cells(::Grid)::IndexCellVector{Int} = @abstractmethod

# @fverdugo Return the encoded extrusion instead?
# If we restrict to polytopes that can be build form
# extrusion tuples, then we cannot accommodate polygonal elements, etc.
# This is a strong limitation IMHO that has to be worked out
"""
Returns the tuple uniquely identifying the Polytope of each cell
"""
function celltypes(::Grid{D,Z})::IndexCellValue{NTuple{Z}} where {D,Z}
  @abstractmethod
end

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

end # module Geometry
