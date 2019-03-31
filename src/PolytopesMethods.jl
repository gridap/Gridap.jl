"""
Constructs a `Polytope` given the type of extrusion, i.e., 1 for "hex" extrusion
and 0 for "tet" extrusion
"""
function Polytope(extrusion::PointInt{D}) where D
  zerop = PointInt{D}(zeros(Int64,D))
	pol_nfs_dim = polytopenfaces(zerop, extrusion)
	pol_nfs = pol_nfs_dim[1]
	pol_dim = pol_nfs_dim[2]
	nfs_id = Dict(nf => i for (i,nf) in enumerate(pol_nfs))
	num_nfs = length(nfs_id)
	nf_nfs_dim = polytopemesh(pol_nfs, nfs_id)
	nf_nfs = nf_nfs_dim[1]; nf_dim = nf_nfs_dim[2]
	Polytope{D}(extrusion, pol_nfs, nf_nfs, nf_dim)
end

"""
Provides the number of n-faces of a polytope
"""
numnftypes(polytope::Polytope) = 2^dim(polytope)


function polytopemesh(nfaces,nfaceid)
	num_nfs = length(nfaces)
	nfnfs = Vector{Vector{Int64}}(undef,num_nfs)
	nfnfs_dim = Vector{Vector{UnitRange{Int64}}}(undef,num_nfs)
	for (inf,nf) in enumerate(nfaces)
		nfs_dim_nf = polytopenfaces(nf.anchor, nf.extrusion)
		# nfs_dim_nf = nfaceboundary!(nf.anchor, zerop, nf.extrusion, true, nfaces)
    nf_nfs = nfs_dim_nf[1]
		dimnfs = nfs_dim_nf[2]
		nfnfs[inf] = [ get(nfaceid,nf,nf) for nf in nf_nfs ]
		nfnfs_dim[inf] = dimnfs
	end
	return [nfnfs, nfnfs_dim]
end

function nfdim(ext::PointInt{D}) where D
	c= 0
	for i in 1:D
		if (ext[i] > 0)
			c +=1
		end
	end
	return c
end

function polytopenfaces(anchor::PointInt{D}, extrusion::PointInt{D}) where D
	# dnf = sum(extrusion)
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
			# global dim
			# global i
			dim+=1
			dimnfs[dim]=i:iface-1
			i=iface
		end
	end
	dimnfs[dnf+1]=numnfs:numnfs
	# aux=2*ones(Int64,D)
	# offst=[ prod(aux[1:i-1]) for i=1:D]
	# nftype=[offst'*nf_nfs[i].extrusion for i=1:numnfs].+1
	# Polytope{D}(extrusion,nf_nfs,dimnfs,nftype)
	return [nf_nfs,dimnfs]
end

function nfaceboundary!(
	anchor::PointInt{D}, extrusion::PointInt{D}, extend::PointInt{D},
	isanchor::Bool, list) where D
	newext = extend
	list = [list..., NFace{D}(anchor, extrusion)]
	for i = 1:D
		curex = newext[i]
		if (curex > 0) # Perform extension
			func1 = (j -> j==i ? 0 : newext[j])
			newext = PointInt{D}([func1(i) for i=1:D])
			func2 = (j -> j==i ? 1 : 0)
			edim = PointInt{D}([func2(i) for i=1:D])
			if (curex == 1) # Quad extension
				list = nfaceboundary!(anchor+edim, extrusion, newext, false, list)
			elseif (isanchor)
				list = nfaceboundary!(edim, extrusion, newext, false, list)
			end
			list = nfaceboundary!(anchor, extrusion+edim*curex, newext, false, list)
		end
	end
	return list
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
