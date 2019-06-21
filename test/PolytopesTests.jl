##
using Gridap, Test
using Gridap.FieldValues
using Gridap.Polytopes
using Gridap.Polytopes: PointInt

# Developing the change of basis for all n-faces of a polytope
# 1. Given an n-face, determine all rigid-body permutations
# 2. Method that identifies the change of basis required to glue
# together two different orientations
# input: [gn1, ..., gnk] global vertex labels in the base n-face order
# input: [gn1', ..., gnk'] in the cell-based order of the n-face
# output: required change of basis
# Case 1: edges in 2D
# gid = [10,11]
# l2g = [10, 11]
# acceptable permutations = [ [1,2], [2,1] ]
# my_perm = 1
# gid = [11, 10]
# my_perm = 2


# Checking all topologies
##
D=3
# Cube
polytope = Polytope(1,1,1)
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

# Checking _dimfrom_fs_dimto_fs
D = 3
p = Polytope(1,1,1)
nf_vs = Gridap.Polytopes._dimfrom_fs_dimto_fs(p,2,0)
@test length(nf_vs) == 6
@test nf_vs[end] == [2,4,6,8]
p = Polytope(1,2,2)
nf_vs = Gridap.Polytopes._dimfrom_fs_dimto_fs(p,2,0)
@test length(nf_vs) == 4
@test nf_vs[end] == [2,3,4]


# Edge
##
D = 3
p = Polytope(1,1,1)
anc = p.nfaces[end].anchor
ext = p.nfaces[end].extrusion

p.nfaces
# vertices coordinates
vs_r = p.nf_dim[end][1]
vs = p.nfaces[vs_r]
num_vs = length(vs)
vs_cs = Vector{Vector{Int}}(undef, num_vs)
for (i_v,v) in enumerate(vs)
  vs_cs[i_v] = v.anchor
end
perm = zeros(Int64, D*num_vs, D)
c = 1
for (i_v,v) in enumerate(vs)
  global(c)
  v_cs = vs_cs[i_v]
  for idim = 1:D
    perm[c, 1] = i_v
    auxv = copy(v_cs)
    v_cs[idim] == 0 ? auxv[idim] = 1 : auxv[idim] = 0
    for k = 1:num_vs
      if vs_cs[k] == auxv
        perm[c, idim+1] = k
        c += 1
        exit
      end
    end
  end
end
perm
##



1 2
2 1

1 -> 2 first axis
1 -> 3 second axis

if placed in 2
2 -> 1 ok
2 -> 3 not possible
2 -> 4 ok

we must do it with i,j !

1 = (0,0)
2 = (1,0)
3 = (0,1)
4 = (1,1)

Thus, at 2
we can choose the first axis to be e1 or e2
If e1, since x2(1) = 1, +(-1,0) -> (0,0)
Next e2, since x2(2) = 0, +(0,1) -> (1,1)
It provides 2,1 and 2,4 all what is needed to determine orientations
t


1 2 3 4
1 3 2 4
2 4 1 3
2 1 4 3
3 4 1 2
3 1 4 2
4 3 2 1
4 2 3 1


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
