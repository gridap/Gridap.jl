module Meshes

using Gridap # @santiagobadia : to be eliminated
using Gridap.Polytopes
using Base.Cartesian
using Gridap.CellValues
using Gridap.CellValues.Wrappers
using Gridap.FieldValues
using Gridap.Helpers

export gidscellxtype, gidscellxtypefast
export flip
export LexIndexSet
export lex2int
export cartesianindexmatrixoffset!
export Mesh, StructHexMesh


# Abstract types and interfaces

abstract type Mesh{D} end

cellvefs(::Mesh)::CellVector{Int} = @abstractmethod

vefcells(::Mesh)::CellVector{Int} = @abstractmethod

# Concrete structs

struct StructHexMesh{D} <: Mesh{D}
  # Here I am computing the cell vefs, and store the result in a big array,
  # probably better to compute it on the fly on demand. Analogously for the
  # dual mesh, for any relation of nface dim against nface dim. For the moment,
  # I am creating a mesh using the same structures as unstructured meshes.
  # A very light weight structured mesh in which arrays are replaced by
  # functions, probably using an AbstractArray would be possible.
  polytope::Polytope{D}
  cellvefs::CellVector{Int}
  vefcells::CellVector{Int}
  # coordinates::Array{Point{D}}
end

cellvefs(this::StructHexMesh)::CellVector{Int} = this.cellvefs

vefcells(this::StructHexMesh)::CellVector{Int} = this.vefcells

struct LexIndexSet
    range::Vector{Int64}
    offset::Vector{Int64}
    function LexIndexSet(range::Vector{Int64})
        offset = [ prod(range[1:i-1]) for i=1:length(range)]
    return new(range,offset)
  end
end

# Methods

include("MeshesMethods.jl")

end # module Meshes
