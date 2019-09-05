module PolytopesTests

##
using Gridap, Test

##
# Eliminating the NodeArray struct
##
# Checking Polytope and NFace APIs
D = 3

p = Polytope(1,1,1)

nf = nfaces(p)[9]

@test dim(p) == 3

@test extrusion(p).array == [1,1,1]

@test length(nfaces(p)) == 27

@test nf_nfs(p)[end] == [i for i in 1:27]

@test nf_dim(p)[end][1] == 1:8

@test dim(nf) == 1

@test space_dim(nf) == 3

@test anchor(nf).array == [0,0,0]

@test vertices_coordinates(p)[8].array == [1.0,1.0,1.0]

@test num_nfaces(p,1) == 12

@test extrusion(nface_ref_polytopes(p)[end]) == extrusion(p)

@test generate_admissible_permutations(p)[end] == [8,7,6,4]

@test length(equidistant_interior_nodes_coordinates(p,4)) == 27

@test equidistant_interior_nodes_coordinates(p,[2,2,2])[1].array == [0.5, 0.5, 0.5]

@test vertices_coordinates(p)[8].array == [1.0,1.0,1.0]

@test facet_normals(p)[1][1].array ≈ [0.0,0.0,-1.0]

## Checking anisotropic order for n-cubes

D = 3

p = Polytope(1,1,1)

order = [2,4,3]

nodes = equidistant_interior_nodes_coordinates(p, order)

@test nodes[3].array ≈ [0.5,0.75,1.0/3.0]

##
# Adding outwards normals
D = 3
p = Polytope(1,1,1)
ns, f_os = Gridap.Polytopes.facet_normals(p)
@test ns[1].array ≈ [0, 0 ,-1]
p = Polytope(1,2,2)
ns, f_os = Gridap.Polytopes.facet_normals(p)
@test ns[1].array ≈ [0, 0 ,-1]
##
# Adding outwards normals
D = 2
p = Polytope(1,1)
ns, f_os = Gridap.Polytopes.facet_normals(p)
@test ns[1].array ≈ [0,-1]
p = Polytope(1,2)
ns, f_os = Gridap.Polytopes.facet_normals(p)
@test ns[1].array ≈ [0,-1]
##
# Adding outwards normals
D = 4
p = Polytope(1,1,1,1)
ns, f_os = Gridap.Polytopes.facet_normals(p)
@test ns[1].array ≈ [0,0,0,-1]
p = Polytope(1,2,2,2)
ns, f_os = Gridap.Polytopes.facet_normals(p)
@test ns[1].array ≈ [0,0,0,-1]
##
D = 3
p = Polytope(1,1,1)
nf_vs = Gridap.Polytopes._dimfrom_fs_dimto_fs(p,2,0)
vs = Gridap.Polytopes.vertices_coordinates(p)
@test length(nf_vs) == 6
@test nf_vs[end] == [2,4,6,8]
p = Polytope(Point(1,2,2))
nf_vs = Gridap.Polytopes._dimfrom_fs_dimto_fs(p,2,0)
@test length(nf_vs) == 4
@test nf_vs[end] == [2,3,4]
##



# using Gridap.Polytopes
p = Polytope(HEX_AXIS,HEX_AXIS)
for dim = 1:4
  p = Polytope(ones(Int,dim)...)
  _order = ones(Int,dim)
  for i in 1:5
    order = Tuple(_order*i)
    # nf = p.nfaces[end]
    vs = Gridap.Polytopes._interior_nodes_int_coords(p,order)
    @test length(vs) == max(0,(i-1)^dim)
  end
  vcs = Gridap.Polytopes.vertices_coordinates(p)
  @test length(vcs) == 2^dim
end
##
a = zeros(Int,4,4)
a[1,2] = 1
a[1,3] = 2
a[1,4] = 3
a[2,3] = 1
a[2,4] = 3
a[3,4] = 1
##
for dim = 1:4
  p = Polytope(2*ones(Int,dim)...)
  _order = ones(Int,dim)
  for i in 1:4
    order = Tuple(_order*i)
    # nf = p.nfaces[end]
    vs = Gridap.Polytopes._interior_nodes_int_coords(p,order)
    @test length(vs) == a[dim,i]
  end
  vcs = Gridap.Polytopes.vertices_coordinates(p)
  @test length(vcs) == dim+1
end
##
p = Polytope(HEX_AXIS, HEX_AXIS, TET_AXIS)
_order = ones(Int,3)
order = Tuple(_order*3)
# nf = p.nfaces[end]
vs = Gridap.Polytopes._interior_nodes_int_coords(p,order)
@test length(vs) == 1
vcs = Gridap.Polytopes.vertices_coordinates(p)
@test length(vcs) == 5
##
pl = Polytope(HEX_AXIS, TET_AXIS, HEX_AXIS)
_order = ones(Int,3)
order = Tuple(_order*3)
# p = pl.nfaces[end]
vs = Gridap.Polytopes._interior_nodes_int_coords(pl,order)
length(vs)
@test length(vs) == 2
vcs = Gridap.Polytopes.vertices_coordinates(pl)
@test length(vcs) == 6
##


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
ext = (1,1,1)
polytope = Polytope(ext)
polytope.extrusion.array.data
x = Vector{Float64}(undef,D)
x .= polytope.nfaces[1].anchor.array
@test length(polytope.nfaces) == 27
@test num_nfaces(polytope,0) == 8
@test num_nfaces(polytope,1) == 12
@test num_nfaces(polytope,2) == 6
@test num_nfaces(polytope,3) == 1
# Pyramid
ext = (1,1,2)
polytope = Polytope(ext)
@test length(polytope.nfaces) == 19
# Prysm
ext = (1,2,1)
polytope = Polytope(ext)
@test length(polytope.nfaces) == 21
ext = Point(1,2,2)
polytope = Polytope(ext)
@test length(polytope.nfaces) == 15
##

# Test to check nodes on the closure of an nface of a polytope
##
D=3
orders=[2,3,2]
ext = Point(1,1,1)
polytope = Polytope(ext)
nodes, nfacenodes = Gridap.RefFEs._high_order_lagrangian_nodes_polytope(polytope,orders)
cv1 = Gridap.CellValuesGallery.CellValueFromArray(polytope.nf_nfs) # cell to index
cv2 = Gridap.CellValuesGallery.CellValueFromArray(nfacenodes)
closurenfacenodes =  Gridap.CellValuesGallery.CellVectorByComposition(cv1,cv2)
@test length(closurenfacenodes[end-1])==12
@test closurenfacenodes[end-1][end-1]==33
##

# Similar test on the (open for dim > 0) nface
##
D=3
orders=[2,3,2]
ext = Point(1,1,1)
polytope = Polytope(ext)
nodes, nfacenodes = Gridap.RefFEs._high_order_lagrangian_nodes_polytope(polytope,orders)
nfacenodes
@test length(nfacenodes[end-1])==2
@test nfacenodes[end-1][end-1]==33
##

# Test to check the node coordinates
##
D=3
orders=[2,3,4]
ext = Point(1,1,1)
polytope = Polytope(ext)
nodes, nfacenodes = Gridap.RefFEs._high_order_lagrangian_nodes_polytope(polytope,orders)
@test length(nodes)==60
@test nodes[33] ≈ Point(0.5, 1.0/3.0, 0.0)
cv1 = Gridap.CellValuesGallery.CellValueFromArray(polytope.nf_nfs) # cell to index
cv2 = Gridap.CellValuesGallery.CellValueFromArray(nfacenodes)
closurenfacenodes =  Gridap.CellValuesGallery.CellVectorByComposition(cv1,cv2)
nfacenodes = closurenfacenodes[end-1]
coords = nodes[nfacenodes,:]
fco = i -> coords[i][1]
@test (prod(fco(i) for i=1:length(coords))==1)
##

# Checking _dimfrom_fs_dimto_fs
D = 3
p = Polytope(1,1,1)
nf_vs = Gridap.Polytopes._dimfrom_fs_dimto_fs(p,2,0)
@test length(nf_vs) == 6
@test nf_vs[end] == [2,4,6,8]
p = Polytope(Point(1,2,2))
nf_vs = Gridap.Polytopes._dimfrom_fs_dimto_fs(p,2,0)
@test length(nf_vs) == 4
@test nf_vs[end] == [2,3,4]

# Edge
##
p = Polytope(1,1,2)
perm_p = Gridap.Polytopes.generate_admissible_permutations(p)
@test length(perm_p) == 48
# @show perm_p
p = Polytope(1,2,2)
perm_p = Gridap.Polytopes.generate_admissible_permutations(p)
@test length(perm_p) == 24
##

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
# D = 3
# p = Polytope(1,1,1)
# anc = p.nfaces[end].anchor
# ext = p.nfaces[end].ext
#
# p.nfaces
# # vertices coordinates
# vs_r = p.nf_dim[end][1]
# vs = p.nfaces[vs_r]
# num_vs = length(vs)
# vs_cs = Vector{Vector{Int}}(undef, num_vs)
# for (i_v,v) in enumerate(vs)
#   vs_cs[i_v] = v.anchor
# end
# perm = zeros(Int64, D*num_vs, D)
# c = 1
# for (i_v,v) in enumerate(vs)
#   global(c)
#   v_cs = vs_cs[i_v]
#   for idim = 1:D
#     perm[c, 1] = i_v
#     auxv = copy(v_cs)
#     v_cs[idim] == 0 ? auxv[idim] = 1 : auxv[idim] = 0
#     for k = 1:num_vs
#       if vs_cs[k] == auxv
#         perm[c, idim+1] = k
#         c += 1
#         exit
#       end
#     end
#   end
# end
# perm
# ##
#
#
#
# 1 2
# 2 1
#
# 1 -> 2 first axis
# 1 -> 3 second axis
#
# if placed in 2
# 2 -> 1 ok
# 2 -> 3 not possible
# 2 -> 4 ok
#
# we must do it with i,j !
#
# 1 = (0,0)
# 2 = (1,0)
# 3 = (0,1)
# 4 = (1,1)
#
# Thus, at 2
# we can choose the first axis to be e1 or e2
# If e1, since x2(1) = 1, +(-1,0) -> (0,0)
# Next e2, since x2(2) = 0, +(0,1) -> (1,1)
# It provides 2,1 and 2,4 all what is needed to determine orientations
# t
#
#
# 1 2 3 4
# 1 3 2 4
# 2 4 1 3
# 2 1 4 3
# 3 4 1 2
# 3 1 4 2
# 4 3 2 1
# 4 2 3 1
#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# Test to check the views of n-face set for a given n
##
# D=3
# ext = Point{D,Int}(1,1,1)
# polytope = Polytope(ext)
# for j=1:length(polytope.ext)+1
#   for i=1:length(polytope.nfaces[polytope.dimnfs[j]])
#     @test (sum(polytope.nfaces[polytope.dimnfs[j]][i].ext)==j-1)
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


end # module
