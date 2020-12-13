module MultiFieldFESpacesTests

using BlockArrays
using FillArrays
using Gridap.Arrays
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
using Gridap.Integration
using Gridap.ReferenceFEs
using Gridap.CellData
using Test

using Gridap.MultiField
using Gridap.Arrays: BlockArrayCooMap

order = 2

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

trian = get_triangulation(model)
degree = order
quad = CellQuadrature(trian,degree)

V = TestFESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:H1)
Q = TestFESpace(model,ReferenceFE(lagrangian,Float64,order-1),conformity=:L2)

U = TrialFESpace(V)
P = TrialFESpace(Q)

multi_field_style = ConsecutiveMultiFieldStyle()

Y = MultiFieldFESpace(Vector{Float64},[V,Q],multi_field_style)
X = MultiFieldFESpace(Vector{Float64},[U,P],multi_field_style)

@test num_free_dofs(X) == num_free_dofs(U) + num_free_dofs(P)
@test num_free_dofs(X) == num_free_dofs(Y)

dy = get_cell_shapefuns(Y)
dv, dq = dy

dx = get_cell_shapefuns_trial(X)
du, dp = dx

cellmat = integrate(dv*du,quad)
cellvec = integrate(dv*2,quad)
cellids = get_cell_to_bgcell(trian)
cellmatvec = pair_arrays(cellmat,cellvec)
@test isa(cellmat, LazyArray{<:Fill{<:BlockArrayCooMap}})
@test is_nonzero_block(cellmat[1],1,1)
@test is_zero_block(cellmat[1],1,2)
@test isa(cellvec, LazyArray{<:Fill{<:BlockArrayCooMap}})
@test is_nonzero_block(cellvec[1],1)
@test is_zero_block(cellvec[1],2)

matvecdata = (cellmatvec,cellids,cellids)
matdata = (cellmat,cellids,cellids)
vecdata = (cellvec,cellids)

free_values = rand(num_free_dofs(X))
xh = FEFunction(X,free_values)
test_fe_function(xh)
@test isa(xh,FEFunction)
uh, ph = xh
@test isa(uh,FEFunction)
@test isa(ph,FEFunction)

cell_isconstr = get_cell_isconstrained(X)
@test cell_isconstr == Fill(false,num_cells(model))

cell_constr = get_cell_constraints(X)
@test isa(cell_constr,LazyArray{<:Fill{<:BlockArrayCooMap}})

cell_dof_ids = get_cell_dof_ids(X)
@test isa(cell_dof_ids,LazyArray{<:Fill{<:BlockArrayCooMap}})

cf = CellField(X,get_cell_dof_ids(X))
@test isa(cf,MultiFieldCellField)

test_fe_space(X,matvecdata,matdata,vecdata)
test_fe_space(Y,matvecdata,matdata,vecdata)

#using Gridap.Visualization
#writevtk(trian,"trian";nsubcells=30,cellfields=["uh" => uh, "ph"=> ph])

f(x) = sin(4*pi*(x[1]-x[2]^2))+1
fh = interpolate([f,f],X)
fh = interpolate_everywhere([f,f],X)
fh = interpolate_dirichlet([f,f],X)

end # module
