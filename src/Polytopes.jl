module Polytopes

using StaticArrays
using Base.Cartesian
using Numa #@fverdugo to be eliminated
using Numa.FieldValues

export Polytope
export NodesArray
export NFace
export dim, numnftypes

# Abstract types and interfaces

# Concrete structs

const PointInt{D} = SVector{D,Int64} where D
# @santiagobadia : Probably add Type of coordinates in Point{D}
# @santiagobadia : I will re-think the NodeArray when I have at my disposal
# the geomap on n-faces, etc. And a clearer definition of the mesh object
# to discuss with @fverdugo

"""
n-face of the polytope, i.e., any polytope of lower dimension `N` representing
its boundary and the polytope itself (for `N` equal to the space dimension `D`)
"""
struct NFace{D}
  anchor::PointInt{D}
  extrusion::PointInt{D}
end

"""
Aggregation of all n-faces that compose the polytope boundary and the polytope
itself, the classification of n-faces with respect to their dimension and type
"""
struct Polytope{D}
  extrusion::PointInt{D}
  nfaces::Vector{NFace}
	nf_nfs::Vector{Vector{Int64}}
	nf_dim::Vector{Vector{UnitRange{Int64}}}
end

"""
Array of nodes for a given polytope and order
"""
struct NodesArray{D}
  coordinates::Vector{Point{D}}
  nfacenodes::Array{Array{Int64,1},1}
  closurenfacenodes::Array{Array{Int64,1},1}
  # @santiagobadia : To be changed to points
end

# Methods

include("PolytopesMethods.jl")

end # module Polytopes
