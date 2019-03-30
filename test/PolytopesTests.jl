##
using Numa, Test
using Numa.FieldValues
using Numa.Polytopes
using Numa.Polytopes: PointInt
using Numa.Polytopes: createnfaces!
##

# Designing new polytope
##
D = 3
extrusion = PointInt{D}(1,1,1)
polytope = Polytope(extrusion)
nfaceid = Dict(nf => i for (i,nf) in enumerate(polytope.nfaces))
numnfs = length(nfaceid)
i=numnfs
nfnfs = Vector{Vector{Int64}}(undef,numnfs)
nfnfs_dim = Vector{Vector{UnitRange{Int64}}}(undef,numnfs)
nf=polytope.nfaces[i]
dnf = sum(nf.extrusion)
nfaces=[]
nf.anchor
nf.extrusion
zerop = PointInt{D}(zeros(Int64,D))
##
nfaces = []
topology = nf.extrusion*2
nfaces = createnfaces!(nf.anchor, zerop, nf.extrusion, nfaces)
nfaces2 = []
nfaces2 = createnfacesnew!(nf.anchor, zerop, nf.extrusion, true, topology, nfaces2)
nfaces == nfaces2
[sort!(nfaces, by = x -> x.anchor[i]) for i=1:length(extrusion)]
[sort!(nfaces, by = x -> x.extrusion[i]) for i=1:length(extrusion)]
[sort!(nfaces, by = x -> sum(x.extrusion))]
numnfs = length(nfaces)
nfsdim = [sum(nfaces[i].extrusion) for i=1:numnfs]
dimnfs = Array{UnitRange{Int64},1}(undef,dnf+1)
dim=0; i=1
for iface=1:numnfs
	if (nfsdim[iface]>dim)
		global dim
		global i
		dim+=1
		dimnfs[dim]=i:iface-1
		i=iface
	end
end
dimnfs[dnf+1]=numnfs:numnfs
nfnfs[i] = [ get(nfaceid,nf,nf) for nf in nfaces ]
nfnfs_dim[i] = dimnfs

function createnfacesnew!(
	anchor::PointInt{D}, extrusion::PointInt{D}, extend::PointInt{D},
	isanchor::Bool, topology::PointInt{D}, list) where D
	newext = copy(extend) # extend cannot be modified for recursion
	a = NFace{D}(anchor, extrusion)
	list = [list...,a]
	for i = 1:D
		if (newext[i] == 1) # Perform extension
			func = (j -> j==i ? 0 : newext[j])
			newext = PointInt{D}([func(i) for i=1:D])
			func = (j -> j==i ? 1 : 0)
			edim = PointInt{D}([func(i) for i=1:D])
			if (topology[i] == 1) # Quad extension
				# newext2 = PointInt{D}([func(i) for i=1:D])
				list = createnfacesnew!(anchor+edim, extrusion, newext, false, topology, list)
			elseif (isanchor)
				func = (j -> j==i ? 0 : newext[j])
				newext = PointInt{D}([func(i) for i=1:D])
				list = createnfacesnew!(edim, extrusion, newext, false, topology, list)
			end
			list = createnfacesnew!(anchor, extrusion+edim, newext, false, topology, list)
		end
	end
	# Create +e_j only once !!! @santiagobadia
	return list
end






##

for (i,nf) in enumerate(polytope.nfaces)
	nfaces = []
	createnfaces!(nf.anchor, nf.extrusion, nf.extrusion, nfaces)
	println("nface ",i)
	println("composition ",nfaces)
end

##

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
