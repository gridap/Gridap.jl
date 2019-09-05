
# 1) Clean constructors without D or T

# Create dofbasis using node array for Lagrangian FEs

# Create BasisWithChangeOfBasis
# i.e., CanonicalBasis given DOFs

# nfacetoowndofs

# D = 1
#


D = 2
orders=[1,1]
orders=[2,3]
p = Polytope(1,1)
reffe = LagrangianRefFE{2,Float64}(p,orders)
#
p = Polytope(1,2)
# end
