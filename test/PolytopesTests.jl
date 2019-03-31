##
using Numa, Test
using Numa.FieldValues
using Numa.Polytopes
using Numa.Polytopes: PointInt
using Numa.Polytopes: createnfaces!
##



struct NewPolytope{D}
  extrusion::PointInt{D}
  nfaces::Vector{NFace}
	nf_nfs::Vector{Vector{Int64}}
	nf_dim::Vector{Vector{UnitRange{Int64}}}
end

function NewPolytope(extrusion::PointInt{D}) where D
  zerop = PointInt{D}(zeros(Int64,D))
	pol_nfs_dim = polytopenfaces(zerop, extrusion)
	pol_nfs = pol_nfs_dim[1]
	pol_dim = pol_nfs_dim[2]
	nfs_id = Dict(nf => i for (i,nf) in enumerate(pol_nfs))
	num_nfs = length(nfs_id)
	nf_nfs_dim = polytopemeshnew(pol_nfs, nfs_id)
	nf_nfs = nf_nfs_dim[1]; nf_dim = nf_nfs_dim[2]
	NewPolytope{D}(extrusion, pol_nfs, nf_nfs, nf_dim)
end

function polytopemeshnew(nfaces,nfaceid)
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

function polytopenfaces(anchor, extrusion::PointInt{D})
	dnf = sum(extrusion)
	zerop = PointInt{D}(zeros(Int64,D))
	nf_nfs = []
	nf_nfs = nfaceboundary!(anchor, zerop, extrusion, true, nf_nfs)
	[sort!(nf_nfs, by = x -> x.anchor[i]) for i=1:length(extrusion)]
	[sort!(nf_nfs, by = x -> x.extrusion[i]) for i=1:length(extrusion)]
	[sort!(nf_nfs, by = x -> sum(x.extrusion))]
	numnfs = length(nf_nfs)
	nfsdim = [sum(nf_nfs[i].extrusion) for i=1:numnfs]
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
			list = nfaceboundary!(anchor, extrusion+edim, newext, false, list)
		end
	end
	return list
end
##
# function polytopemesh!(nfaces,nfaceid,nfnfs,nfnfs_dim)
# 	for (inf,nf) in enumerate(nfaces)
# 		dnf = sum(nf.extrusion)
# 		zerop = PointInt{D}(zeros(Int64,D))
# 		nf_nfs = []
# 		nf_nfs = nfaceboundary!(nf.anchor, zerop, nf.extrusion, true, nf_nfs)
# 		[sort!(nf_nfs, by = x -> x.anchor[i]) for i=1:length(extrusion)]
# 		[sort!(nf_nfs, by = x -> x.extrusion[i]) for i=1:length(extrusion)]
# 		[sort!(nf_nfs, by = x -> sum(x.extrusion))]
# 		numnfs = length(nf_nfs)
# 		nfsdim = [sum(nf_nfs[i].extrusion) for i=1:numnfs]
# 		dimnfs = Array{UnitRange{Int64},1}(undef,dnf+1)
# 		dim=0; i=1
# 		for iface=1:numnfs
# 			if (nfsdim[iface]>dim)
# 				# global dim
# 				# global i
# 				dim+=1
# 				dimnfs[dim]=i:iface-1
# 				i=iface
# 			end
# 		end
# 		dimnfs[dnf+1]=numnfs:numnfs
# 		nfnfs[inf] = [ get(nfaceid,nf,nf) for nf in nf_nfs ]
# 		nfnfs_dim[inf] = dimnfs
# 	end
# end




# Designing new polytope
##
D = 3
extrusion = PointInt{D}(1,1,1)
polytope = Polytope(extrusion)
##
# function Polytope(extrusion::PointInt{D}) where D
zerop = PointInt{D}(zeros(Int64,D))
pol_nfs_dim = polytopenfaces(zerop, extrusion)
pol_nfs = pol_nfs_dim[1]
pol_dim = pol_nfs_dim[2]
nfs_id = Dict(nf => i for (i,nf) in enumerate(pol_nfs))
num_nfs = length(nfs_id)
##
nf_nfs_dim = polytopemeshnew(pol_nfs, nfs_id)
nf_nfs = nf_nfs_dim[1]; nf_dim = nf_nfs_dim[2]
##
new_poly = NewPolytope(extrusion)


# Test to check nodes on the closure of an nface of a polytope
##
D=3
orders=[2,3,2]
extrusion = PointInt{D}(1,1,1)
polytope = Polytope(extrusion)
nodes = NodesArray(polytope,orders)
@test length(nodes.closurenfacenodes[end-1])==12
@test nodes.closurenfacenodes[end-1][end-1]==33
##

# Similar test on the (open for dim > 0) nface
##
D=3
orders=[2,3,2]
extrusion = PointInt{D}(1,1,1)
polytope = Polytope(extrusion)
nodes = NodesArray(polytope,orders)
@test length(nodes.nfacenodes[end-1])==2
@test nodes.nfacenodes[end-1][end-1]==18
##

# Test to check the node coordinates
##
D=3
orders=[2,3,4]
extrusion = PointInt{D}(1,1,1)
polytope = Polytope(extrusion)
nodes = NodesArray(polytope,orders)
@test length(nodes.coordinates)==60
@test nodes.coordinates[33] ≈ [1.0, 1.0/3.0, 0.0]
nfacenodes = nodes.closurenfacenodes[end-1]
coords = nodes.coordinates[nfacenodes,:]
f = i -> coords[i][1]
@test (prod(f(i) for i=1:length(coords))==1)
##

# Test to check the views of n-face set for a given n
##
D=3
extrusion = PointInt{D}(1,1,1)
polytope = Polytope(extrusion)
for j=1:length(polytope.extrusion)+1
	for i=1:length(polytope.nfaces[polytope.dimnfs[j]])
		@test (sum(polytope.nfaces[polytope.dimnfs[j]][i].extrusion)==j-1)
	end
end
##
