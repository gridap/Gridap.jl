module Meshes

using Numa # @santiagobadia : to be eliminated
using Numa.Polytopes
using Base.Cartesian

export gidscellxtype, gidscellxtypefast
export flip
export LexIndexSet
export lex2int
export cartesianindexmatrixoffset!
export Mesh, StructHexMesh


# Abstract types and interfaces

abstract type Mesh end

# Concrete structs

struct StructHexMesh <: Mesh
	# Here I am computing the cel vefs, and store the result in a big array, probably better to compute it on the fly on demand. Analogously for the dual mesh, for any relation of nface dim against nface dim. For the moment, I am creating a mesh using the same structures as unstructured meshes. A very light weight structured mesh in which arrays are replaced by functions, probably using an AbstractArray would be possible.
	cellvefs::Array{Array{Int64,1},1}
	vefcells::Array{Array{Int64,1},1}
	coordinates::Array{Float64,2}
	polytope::Polytope
end

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
