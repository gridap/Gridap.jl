module DOFBasesTests

using Test
using Gridap

fun(x) = x[1]+x[2]

D = 2
T = Float64
field = AnalyticalField(fun,D)

basis = MonomialBasis(T,(1,1))

polytope = Polytope(HEX_AXIS,HEX_AXIS)
nodes, nfacenodes = Gridap.RefFEs._high_order_lagrangian_nodes_polytope(polytope,[1,1])
nodes
dofbasis = LagrangianDOFBasis(T,nodes)

field_dofs = [0.0, 1.0, 1.0, 2.0]

basis_dofs = [
 1.0 1.0 1.0 1.0;
 0.0 1.0 0.0 1.0;
 0.0 0.0 1.0 1.0;
 0.0 0.0 0.0 1.0]

test_dof_basis(dofbasis,field,basis,field_dofs,basis_dofs)

T = VectorValue{2,Float64}

fun(x) = VectorValue(x[1]+x[2],x[1])

field = AnalyticalField(fun,D)
basis = MonomialBasis(T,(1,1))

nodes, nfacenodes = Gridap.RefFEs._high_order_lagrangian_nodes_polytope(polytope,[1,1])
dofbasis = LagrangianDOFBasis(T,nodes)

# # Component major
#field_dofs = [0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 2.0, 1.0]

# Node major
field_dofs = [0.0, 1.0, 1.0, 2.0, 0.0, 1.0, 0.0, 1.0]

# Node major
basis_dofs = [
  1.0 1.0 1.0 1.0 0.0 0.0 0.0 0.0;
  0.0 0.0 0.0 0.0 1.0 1.0 1.0 1.0;
  0.0 1.0 0.0 1.0 0.0 0.0 0.0 0.0;
  0.0 0.0 0.0 0.0 0.0 1.0 0.0 1.0;
  0.0 0.0 1.0 1.0 0.0 0.0 0.0 0.0;
  0.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0;
  0.0 0.0 0.0 1.0 0.0 0.0 0.0 0.0;
  0.0 0.0 0.0 0.0 0.0 0.0 0.0 1.0]

test_dof_basis(dofbasis,field,basis,field_dofs,basis_dofs)

@test dofbasis.dof_to_node == [1, 2, 3, 4, 1, 2, 3, 4]
@test dofbasis.dof_to_comp == [1, 1, 1, 1, 2, 2, 2, 2]
@test dofbasis.node_and_comp_to_dof == [1 5; 2 6; 3 7; 4 8]

end # module
