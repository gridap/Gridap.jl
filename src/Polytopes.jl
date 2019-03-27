export Polytope
export NodesArray
export NFace
export dim, numnftypes

const PointInt{D} = SVector{D,Int64} where D
# @santiagobadia : Probably add Type of coordinates in Point{D}

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
  dimnfs::Vector{UnitRange{Int64} }
  nftype::Vector{Int64}
end

function Polytope(extrusion::PointInt{D}) where D
  zerop = PointInt{D}(zeros(Int64,D))
  nfaces = []
  nfaces = createnfaces!(zerop, zerop, extrusion, nfaces)
  [sort!(nfaces, by = x -> x.anchor[i]) for i=1:length(extrusion)]
  [sort!(nfaces, by = x -> x.extrusion[i]) for i=1:length(extrusion)]
  [sort!(nfaces, by = x -> sum(x.extrusion))]
  numnfs = length(nfaces)
  nfsdim = [sum(nfaces[i].extrusion) for i=1:numnfs]
  dimnfs = Array{UnitRange{Int64},1}(undef,D+1)
  dim=0; i=1
  for iface=1:numnfs
    if (nfsdim[iface]>dim)
      # global dim
      # global i
      dim+=1
      dimnfs[dim]=i:iface-1
      i=iface
    end
  end
  dimnfs[D+1]=numnfs:numnfs
  aux=2*ones(Int64,D)
  offst=[ prod(aux[1:i-1]) for i=1:D]
  nftype=[offst'*nfaces[i].extrusion for i=1:numnfs].+1
  Polytope{D}(extrusion,nfaces,dimnfs,nftype)
end

# dim(polytope::Polytope) = length(polytope.extrusion)
numnftypes(polytope::Polytope) = 2^dim(polytope)

# Create list of all nfaces of a polytope sorted by dim first
function createnfaces!(
  anchor::PointInt{D},
  extrusion::PointInt{D},
  extend::PointInt{D},
  list) where D
  newext = copy(extend) # extend cannot be modified for recursion
  a = NFace{D}(anchor, extrusion)
  list = [list...,a]
  for i = 1:D
    if (newext[i] > 0)
      func = (j -> j==i ? 1 : 0)
      edim = PointInt{D}([func(i) for i=1:D])
      func = (j -> j==i ? 0 : newext[j])
      newext = PointInt{D}([func(i) for i=1:D])
      list = createnfaces!(anchor+edim, extrusion, newext, list)
      list = createnfaces!(anchor, extrusion+edim, newext, list)
    end
  end
  return list
end

"""
Array of nodes for a given polytope and order
"""
struct NodesArray{D}
  coordinates::Vector{Point{D}}
  nfacenodes::Array{Array{Int64,1},1}
  closurenfacenodes::Array{Array{Int64,1},1}
end

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
  points = Vector{Point{D}}(undef,npoints)
  cid = ntuple(i -> 1:length(coords[i]), D)
  cid = CartesianIndices(cid)
  tpcoor = j -> [ coords[i][j[i]] for i ∈ 1:D]
  for (i,j) ∈ enumerate(cid)
    points[i] = tpcoor(Tuple(j))
  end
  NodesArray(points, nfacenodes, closurenfacenodes)
end

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
