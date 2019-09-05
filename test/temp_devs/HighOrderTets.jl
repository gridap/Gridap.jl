module HighOrderTets

##
using Gridap, Test
# 1) Clean constructors without D or T

# Create dofbasis using node array for Lagrangian FEs

# Create BasisWithChangeOfBasis
# i.e., CanonicalBasis given DOFs

# nfacetoowndofs

# D = 1
#


D = 2
orders=[3,3]
p = Gridap.Polytope(1,1)
reffe = LagrangianRefFE{2,Float64}(p,orders)
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
