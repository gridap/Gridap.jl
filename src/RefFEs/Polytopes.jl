module Polytopes

using Gridap
using StaticArrays
using Base.Cartesian

export Polytope
export NodesArray
export NFace
export HEX_AXIS, TET_AXIS
export num_nfaces

# Module constants

const HEX_AXIS = 1
const TET_AXIS = 2

# Concrete structs and their pubic API

const PointInt{D} = SVector{D,Int64} where D
# @santiagobadia : Probably add Type of coordinates in Point{D} (@fverdugo: Now we have Point{D,T})
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
Constructs a `Polytope` given the type of extrusion, i.e., HEX_AXIS (=1) for "hex" extrusion
and TET_AXIS (=2) for "tet" extrusion
"""
function Polytope(extrusion::PointInt{D}) where D
  zerop = PointInt{D}(zeros(Int64,D))
  pol_nfs_dim = polytopenfaces(zerop, extrusion)
  pol_nfs = pol_nfs_dim[1]; pol_dim = pol_nfs_dim[2]
  nfs_id = Dict(nf => i for (i,nf) in enumerate(pol_nfs))
  nf_nfs_dim = polytopemesh(pol_nfs, nfs_id)
  nf_nfs = nf_nfs_dim[1]; nf_dim = nf_nfs_dim[2]
  Polytope{D}(extrusion, pol_nfs, nf_nfs, nf_dim)
end

function Polytope(extrusion::NTuple{D,Int}) where D
  Polytope(SVector(extrusion...))
end

function Polytope(extrusion::Vararg{Int,D}) where D
  Polytope(extrusion)
end

#@fverdugo add here a public API to access the polytope info instead
# of using the struct fields directly as it is currently done
# in some places of the code.

"""
Returns number of nfaces of dimension dim
"""
function num_nfaces(polytope::Polytope, dim::Integer)
  n = length(polytope.nf_dim)
  k = 0
  for nface in 1:n
    d = length(polytope.nf_dim[nface])-1
    if d == dim
      k +=1
    end
  end
  k
end

"""
Array of nodes for a given polytope and order
"""
struct NodesArray{D}
  coordinates::Vector{Point{D,Float64}}
  nfacenodes::Array{Array{Int64,1},1}
  closurenfacenodes::Array{Array{Int64,1},1}
  # @santiagobadia : To be changed to points
end

"""
Creates an array of nodes `NodesArray` for a given polytope and the order per
dimension
"""
function NodesArray(polytope::Polytope, orders::Array{Int64,1})
  closurenfacenodes = [
  createnodes(polytope.nfaces[i],
  orders) for i=1:length(polytope.nfaces)]
  nfacenodes = [
  createnodes(polytope.nfaces[i],orders,
  isopen=true) for i=1:length(polytope.nfaces)]
  D = length(orders)
  coords = [
  nodescoordinates(orders[i], nodestype="Equispaced") for i=1:D]
  npoints =  prod(ntuple(i -> length(coords[i]), D))
  points = Vector{Point{D,Float64}}(undef,npoints)
  cid = ntuple(i -> 1:length(coords[i]), D)
  cid = CartesianIndices(cid)
  tpcoor = j -> [ coords[i][j[i]] for i ∈ 1:D]
  for (i,j) ∈ enumerate(cid)
    points[i] = tpcoor(Tuple(j))
  end
  NodesArray(points, nfacenodes, closurenfacenodes)
end


# Helpers

# @fverdugo Following the style in julia.Base, I would name with a leading
# underscore all helper functions that are not exposed to the public API

nfdim(a::PointInt{D}) where D = sum([a[i] > 0 ? 1 : 0 for i =1:D ])

"""
Generates the array of n-faces of a polytope
"""
function polytopenfaces(anchor::PointInt{D}, extrusion::PointInt{D}) where D
  dnf = nfdim(extrusion)
  zerop = PointInt{D}(zeros(Int64,D))
  nf_nfs = []
  nf_nfs = nfaceboundary!(anchor, zerop, extrusion, true, nf_nfs)
  [sort!(nf_nfs, by = x -> x.anchor[i]) for i=1:length(extrusion)]
  [sort!(nf_nfs, by = x -> x.extrusion[i]) for i=1:length(extrusion)]
  [sort!(nf_nfs, by = x -> sum(x.extrusion))]
  numnfs = length(nf_nfs)
  nfsdim = [nfdim(nf_nfs[i].extrusion) for i=1:numnfs]
  dimnfs = Array{UnitRange{Int64},1}(undef,dnf+1)
  dim=0; i=1
  for iface=1:numnfs
    if (nfsdim[iface]>dim)
      # global dim; # global i
      dim+=1
      dimnfs[dim]=i:iface-1
      i=iface
    end
  end
  dimnfs[dnf+1]=numnfs:numnfs
  return [nf_nfs,dimnfs]
end

"""
Provides for all n-faces of a polytope the d-faces for 0 <= d <n on its
boundary (e.g., given a face, it provides the set of edges and corners on its
boundary) using the global n-face numbering of the base polytope
"""
function polytopemesh(nfaces::Vector{NFace{D}}, nfaceid::Dict) where D
  num_nfs = length(nfaces)
  nfnfs = Vector{Vector{Int64}}(undef,num_nfs)
  nfnfs_dim = Vector{Vector{UnitRange{Int64}}}(undef,num_nfs)
  for (inf,nf) in enumerate(nfaces)
    nfs_dim_nf = polytopenfaces(nf.anchor, nf.extrusion)
    nf_nfs = nfs_dim_nf[1]; dimnfs = nfs_dim_nf[2]
    nfnfs[inf] = [ get(nfaceid,nf,nf) for nf in nf_nfs ]
    nfnfs_dim[inf] = dimnfs
  end
  return [nfnfs, nfnfs_dim]
end

"""
Generates the list of n-face of a polytope the d-faces for 0 <= d <n on its
boundary
"""
function nfaceboundary!(
  anchor::PointInt{D}, extrusion::PointInt{D},
  extend::PointInt{D}, isanchor::Bool, list) where D
  newext = extend
  list = [list..., NFace{D}(anchor, extrusion)]
  for i = 1:D
    curex = newext[i]
    if (curex > 0) # Perform extension
      func1 = (j -> j==i ? 0 : newext[j])
      newext = PointInt{D}([func1(i) for i=1:D])
      func2 = (j -> j==i ? 1 : 0)
      edim = PointInt{D}([func2(i) for i=1:D])
      func3 = (j -> j >= i ? anchor[j] : 0 )
      tetp = PointInt{D}([func3(i) for i=1:D]) + edim
      if (curex == 1) # Quad extension
        list = nfaceboundary!(anchor+edim, extrusion, newext, false, list)
      elseif (isanchor)
        list = nfaceboundary!(tetp, extrusion, newext, false, list)
      end
      list = nfaceboundary!(anchor, extrusion+edim*curex, newext, false, list)
    end
  end
  return list
end

"""
It provides for every df-face in the polytope all its dt-faces
"""
# @santiagobadia : New method required by @fverdugo
function _dimfrom_fs_dimto_fs(p::Polytope, dim_from::Int, dim_to::Int)
  @assert dim_to <= dim_from
  dim_from +=1; dim_to +=1
  dffs_r = p.nf_dim[end][dim_from]
  dffs_dtfs = Vector{Vector{Int}}(undef, dffs_r[end]-dffs_r[1]+1)
  offs = p.nf_dim[end][dim_to][1]-1
  for (i_dff, dff) in enumerate(dffs_r)
    dff_nfs = p.nf_nfs[dff]
    dff_dtfs_r = p.nf_dim[dff][dim_to]
    dff_dtfs = dff_nfs[dff_dtfs_r]
    dffs_dtfs[i_dff] = dff_dtfs .- offs
    # @santiagobadia : With or without offset ?
  end
  return dffs_dtfs
end

# @santiagobadia : The rest is waiting for a geomap

# Create list of nface nodes with polytope indexing
function createnodes(nface::NFace,orders; isopen::Bool=false)
  spdims=length(nface.extrusion)
  @assert spdims == length(orders) "nface and orders dim must be identical"
  perms=nfaceperms(nface)
  p=perms[1]; pinv=perms[2]
  ordnf=orders[pinv]; ordp1=ordnf.+((isopen) ? d=-1 : d=1)
  A = Array{Array{Int64,1}}(undef,ordp1...)
  cartesianindexmatrix!(A)
  A = hcat(reshape(A, length(A))...)'
  if (isopen) A.+=1 end
  B = Array{Int64,2}(undef,prod(ordp1),spdims)
  for i=1:spdims
    if (p[i]==0) B[:,i].=nface.anchor[i]*orders[i] end
  end
  B[:,pinv] = copy(A)
  offst = orders.+1; offst = [ prod(offst[1:i-1]) for i=1:length(offst)]
  return B*offst.+1
end

# nface to polytope cartesian index change of basis and inverse
function nfaceperms(nface)
  pd = length(nface.extrusion)
  e = nface.extrusion
  c=0; p=[]
  for i=1:pd
    (e[i]!=0) ? (c+=1; p=[p...,c]) : p=[p...,0]
  end
  c=0; pinv=[]
  for i=1:pd
    (p[i]!=0) ? (c+=1; pinv=[pinv...,i]) : 0
  end
  return [p,pinv]
end

# 1-dim node coordinates
function nodescoordinates(order::Int; nodestype::String="Equispaced")
  @assert ((nodestype=="Equispaced") | (nodestype=="Chebyshev")) "Node type not implemented"
  ordp1 = order+1
  if nodestype == "Chebyshev"
    nodescoordinates = [ cos((i-1)*pi/order) for i=1:ordp1]
  elseif nodestype == "Equispaced"
    (order != 0) ? nodescoordinates = (2/order)*[ i-1 for i=1:ordp1].-1 : nodescoordinates = [0.0]
  end
  return nodescoordinates
end

@generated function cartesianindexmatrix!(A::Array{Array{T,M},N}) where {T,M,N}
  quote
    @nloops $N i A begin (@nref $N A i) = [((@ntuple $N i).-1)...] end
  end
end

end # module Polytopes
