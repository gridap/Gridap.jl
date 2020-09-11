module MultiFieldFESpacesTests

using BlockArrays
using FillArrays
using Gridap.Arrays
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
using Gridap.Integration
using Gridap.CellData
using Test

using Gridap.MultiField

order = 2

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

trian = get_triangulation(model)
degree = order
quad = CellQuadrature(trian,degree)
q = get_coordinates(quad)

ref_style = [:reference,:physical]

for ref_st in ref_style

  V = TestFESpace(model=model,order=order,reffe=:Lagrangian,conformity=:H1,valuetype=Float64,dof_space=ref_st)
  Q = TestFESpace(model=model,order=order-1,reffe=:Lagrangian,conformity=:L2,valuetype=Float64,dof_space=ref_st)

  U = TrialFESpace(V)
  P = TrialFESpace(Q)

  multi_field_style = ConsecutiveMultiFieldStyle()

  Y = MultiFieldFESpace([V,Q],multi_field_style)
  X = MultiFieldFESpace([U,P],multi_field_style)

  cell_axes = get_cell_axes(Y)
  @test isa(cell_axes[1][1],BlockedUnitRange)

  cell_axes = get_cell_axes(X)
  @test isa(cell_axes[1][1],BlockedUnitRange)

  @test num_free_dofs(X) == num_free_dofs(U) + num_free_dofs(P)
  @test num_free_dofs(X) == num_free_dofs(Y)

  dy = get_cell_basis(Y)
  @test is_test(dy)
  @test is_a_fe_cell_basis(dy)
  dv, dq = dy
  @test is_a_fe_cell_basis(dv)
  @test is_a_fe_cell_basis(dq)

  dx = get_cell_basis(X)
  @test is_a_fe_cell_basis(dx)
  @test is_trial(dx)
  du, dp = dx
  @test is_a_fe_cell_basis(du)
  @test is_a_fe_cell_basis(dp)

  cellmat = integrate(dv*du,quad)
  cellvec = integrate(dv*2,quad)
  cellids = get_cell_id(trian)
  cellmatvec = pair_arrays(cellmat,cellvec)
  @test isa(cellmat, VectorOfBlockArrayCoo)
  @test is_nonzero_block(cellmat,1,1)
  @test is_zero_block(cellmat,1,2)
  @test isa(cellvec, VectorOfBlockArrayCoo)
  @test is_nonzero_block(cellvec,1)
  @test is_zero_block(cellvec,2)

  matvecdata = (cellmatvec,cellids,cellids)
  matdata = (cellmat,cellids,cellids)
  vecdata = (cellvec,cellids)

  free_values = rand(num_free_dofs(X))
  xh = FEFunction(X,free_values)
  test_fe_function(xh)
  @test is_a_fe_function(xh)
  uh, ph = xh
  @test is_a_fe_function(uh)
  @test is_a_fe_function(ph)

  cell_isconstr = get_cell_isconstrained(X)
  @test cell_isconstr == Fill(false,num_cells(model))

  cell_constr = get_cell_constraints(X)
  @test isa(cell_constr,VectorOfBlockArrayCoo)

  cell_dof_ids = get_cell_dofs(X)
  @test isa(cell_dof_ids,VectorOfBlockArrayCoo)

  cf = CellField(X,get_cell_dofs(X))
  @test isa(cf,MultiFieldCellField)

  test_fe_space(X,matvecdata,matdata,vecdata)
  test_fe_space(Y,matvecdata,matdata,vecdata)

  #using Gridap.Visualization
  #writevtk(trian,"trian";nsubcells=30,cellfields=["uh" => uh, "ph"=> ph])

  f(x) = sin(4*pi*(x[1]-x[2]^2))+1
  fh = interpolate([f,f],X)
  fh = interpolate_everywhere([f,f],X)
  fh = interpolate_dirichlet([f,f],X)
end

end # module
