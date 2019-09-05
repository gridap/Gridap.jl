module HighOrderTets

##
using Gridap, Test
using Gridap.CellValuesGallery
# 1) Clean constructors without D or T

# Create dofbasis using node array for Lagrangian FEs

# Create BasisWithChangeOfBasis
# i.e., CanonicalBasis given DOFs

# nfacetoowndofs

# D = 1
#


# Closure n-face nodes
# Method that given the set of nodes and the nfacedofs, returns the
# nface dofs on the closure of the n-face (only sense for Lagrangian FEs)
##
cell_to_x_l = [2,3,1,3,4,4,3,2,5,4,3,4]
cell_to_x_p = [1,4,4,7,13]
cell_to_x = CellVectorFromDataAndPtrs(cell_to_x_l,cell_to_x_p)
x_to_vals_l = [5,4,1,2,3,6,7,8,9,10]
x_to_vals_p = [1,3,5,6,10,11]
x_to_vals = CellVectorFromDataAndPtrs(x_to_vals_l,x_to_vals_p)
cell_to_vals = CellVectorByComposition(cell_to_x, x_to_vals)
a = [
  [1, 2, 3, 5, 4],
  Int[],
  [3, 6, 7, 8, 9, 6, 7, 8, 9],
  [3, 1, 2, 10, 6, 7, 8, 9, 3, 6, 7, 8, 9]]
##
cell_to_vals
D = 2
orders=[3,3]
p = Gridap.Polytope(1,1)
reffe = LagrangianRefFE{2,Float64}(p,orders)

reffe.nfacedofs
p.nf_nfs
reffe.nfacedofs[p.nf_nfs[i_nf]]

reffe.polytope
cv1 = CellValueFromArray(reffe.polytope.nf_nfs) # cell to index
cv2 = CellValueFromArray(reffe.nfacedofs)
a = CellVectorByComposition(cv1,cv2)




i_nf = 9
nf_nfs = reffe.nfacedofs[p.nf_nfs[i_nf]]


length(orders) == 2
#
p = Polytope(1,2)

orders=[2,3]

nodes, nfacenodes = Gridap.RefFEs._high_order_lagrangian_nodes_polytope(p,orders)

orders = [1,1,1]

if !( all(orders.==1))
    b=10
end
b

end # module
