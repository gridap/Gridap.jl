module CellDofBasesTests

using Test
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Arrays
using Gridap.Fields
using Gridap.CellData

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

order = 1

V = TestFESpace(
  reffe = :Lagrangian,
  conformity = :H1,
  valuetype=Float64,
  order = order,
  model = model,
  dirichlet_tags = [1,6],
  dof_space = :physical)

cell_dof_basis = get_cell_dof_basis(V)
@test is_in_physical_space(cell_dof_basis)

indices = [2,3,3,1]

f = reindex(cell_dof_basis,indices)
@test length(get_array(f)) == length(indices)
@test is_in_physical_space(f)

cell_basis = get_cell_basis(V)
@test is_in_physical_space(cell_basis)

test_cell_dof_basis(cell_dof_basis,cell_basis)

#display(evaluate(cell_dof_basis,cell_basis))

V = TestFESpace(
  reffe = :Lagrangian,
  conformity = :H1,
  valuetype=Float64,
  order = order,
  model = model,
  dirichlet_tags = [1,6])

cell_dof_basis = get_cell_dof_basis(V)
@test is_in_ref_space(cell_dof_basis)

cell_basis = get_cell_basis(V)
@test is_in_ref_space(cell_basis)

test_cell_dof_basis(cell_dof_basis,cell_basis)

#display(evaluate(cell_dof_basis,cell_basis))

end # module
