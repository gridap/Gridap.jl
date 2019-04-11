
# @fverdugo the name of the following types are likely to change

"""
Minimal interface needed to describe a conforming linear discretization
of a domain. E.g., this is the minimum needed to visualize the results.

D is the dimension of the coordinates and Z is the dimension of the cells
"""
abstract type Grid{D,Z} end

function points(::Grid{D})::IndexCellValue{Point{D}} where D
  @abstractmethod
end

cells(::Grid)::IndexCellVector{Int} = @abstractmethod

# @fverdugo Return the encoded extrusion instead?
"""
Returns the tuple uniquely identifying the Polytope of each cell
"""
function celltypes(::Grid{D,Z})::IndexCellValue{NTuple{Z}} where {D,Z}
  @abstractmethod
end


# @fverdugo move this to another place, Polytope.jl?

"""
Encodes the tuple defining a Polytope into an integer
"""
function encode_extrusion(extrusion::NTuple{Z,Int}) where Z
  k = 0
  for (i,v) in enumerate(extrusion)
    k += v*3^i
  end
  k
end

"""
Decodes an integer into a tuple defining a Polytope
"""
function decode_extrusion(i::Int,::Val{Z}) where Z
  @notimplemented
end

