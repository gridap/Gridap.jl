
# @fverdugo the name of the following types are likely to change


#@fverdugo TODO grid could also be high order
"""
D is the dimension of the coordinates and Z is the dimension of the cells
"""
abstract type Grid{D,Z} end # <: Triangulation{Z,D}

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


#@fverdugo put this in another place

struct GridGraphFromData{C<:IndexCellVector{Int},V<:IndexCellVector{Int}} <: GridGraph
  celltovefs::C
  veftocells::V
end

celltovefs(self::GridGraphFromData) = self.celltovefs

veftocells(self::GridGraphFromData) = self.veftocells






