##
using Numa, Test
using Numa.FieldValues
using Numa.Polytopes
using Numa.Polytopes: PointInt
##



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


# @santiagobadia : Add tests for tets

# Designing new polytope
##
# D = 3
# extrusion = PointInt{D}(1,1,1)
# polytope = Polytope(extrusion)
# ##
# # function Polytope(extrusion::PointInt{D}) where D
# zerop = PointInt{D}(zeros(Int64,D))
# pol_nfs_dim = polytopenfaces(zerop, extrusion)
# pol_nfs = pol_nfs_dim[1]
# pol_dim = pol_nfs_dim[2]
# nfs_id = Dict(nf => i for (i,nf) in enumerate(pol_nfs))
# num_nfs = length(nfs_id)
# ##
# nf_nfs_dim = polytopemeshnew(pol_nfs, nfs_id)
# nf_nfs = nf_nfs_dim[1]; nf_dim = nf_nfs_dim[2]
# ##
# new_poly = NewPolytope(extrusion)


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
# D=3
# extrusion = PointInt{D}(1,1,1)
# polytope = Polytope(extrusion)
# for j=1:length(polytope.extrusion)+1
# 	for i=1:length(polytope.nfaces[polytope.dimnfs[j]])
# 		@test (sum(polytope.nfaces[polytope.dimnfs[j]][i].extrusion)==j-1)
# 	end
# end
##
