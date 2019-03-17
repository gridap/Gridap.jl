using Numa, Test

# Test to check nodes on the closure of an nface of a polytope
##
orders=[2,3,2]
polytope = Polytope([1,1,1])
nodes = NodesArray(polytope,orders)
@test length(nodes.closurenfacenodes[end-1])==12
@test nodes.closurenfacenodes[end-1][end-1]==33
##

# Similar test on the (open for dim > 0) nface
##
orders=[2,3,2]
polytope = Polytope([1,1,1])
nodes = NodesArray(polytope,orders)
@test length(nodes.nfacenodes[end-1])==2
@test nodes.nfacenodes[end-1][end-1]==18
##

# Test to check the node coordinates
##
orders=[2,3,4]
polytope = Polytope([1,1,1])
nodes = NodesArray(polytope,orders)
@test length(nodes.coordinates)==180
@test nodes.coordinates[33,:] ≈ [1.0, 1.0/3.0, 0.0]
nfacenodes = nodes.closurenfacenodes[end-1]
coords = nodes.coordinates[nfacenodes,:]
@test (prod(coords[:,1])==1)
##

# Test to check the views of n-face set for a given n
##
polytope = Polytope([1,1,1])
for j=1:length(polytope.extrusion)+1
	for i=1:length(polytope.nfaces[polytope.dimnfs[j]])
		@test (sum(polytope.nfaces[polytope.dimnfs[j]][i].extrusion)==j-1)
	end
end
##
