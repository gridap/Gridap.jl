using Numa, Test
using Numa.Polytopes
using Numa.Polytopes: PointInt

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
