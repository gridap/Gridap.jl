##
using Gridap, Test
using Gridap.FieldValues
using Gridap.Polytopes
using Gridap.Polytopes: PointInt

# Checking all topologies
##
D=3
# Cube
extrusion = PointInt{D}(1,1,1)
polytope = Polytope(extrusion)
@test length(polytope.nfaces) == 27
@test num_nfaces(polytope,0) == 8
@test num_nfaces(polytope,1) == 12
@test num_nfaces(polytope,2) == 6
@test num_nfaces(polytope,3) == 1
# Pyramid
extrusion = PointInt{D}(1,1,2)
polytope = Polytope(extrusion)
@test length(polytope.nfaces) == 19
# Prysm
extrusion = PointInt{D}(1,2,1)
polytope = Polytope(extrusion)
@test length(polytope.nfaces) == 21
extrusion = PointInt{D}(1,2,2)
polytope = Polytope(extrusion)
@test length(polytope.nfaces) == 15
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
fco = i -> coords[i][1]
@test (prod(fco(i) for i=1:length(coords))==1)
##

# Test to check the views of n-face set for a given n
##
# D=3
# extrusion = PointInt{D}(1,1,1)
# polytope = Polytope(extrusion)
# for j=1:length(polytope.extrusion)+1
#   for i=1:length(polytope.nfaces[polytope.dimnfs[j]])
#     @test (sum(polytope.nfaces[polytope.dimnfs[j]][i].extrusion)==j-1)
#   end
# end
##
# using Gridap.Vtkio
# D = 3
# wedge = Polytope(1,2,1)
# @show wedge.nf_nfs
# writevtk(wedge,"wedge")
# pyramid = Polytope(1,1,2)
# writevtk(pyramid,"pyramid")
# hex = Polytope(1,1,1)
# writevtk(hex,"hex")
# tet = Polytope(1,2,2)
# writevtk(tet,"tet")
