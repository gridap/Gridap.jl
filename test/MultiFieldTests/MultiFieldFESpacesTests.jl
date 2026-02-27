module MultiFieldFESpacesTests

using FillArrays
using BlockArrays
using Gridap.Arrays
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Algebra
using Test

using Gridap.MultiField

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

X = MultiFieldFESpace(Vector{Float64},[U,P])

mfs = ConsecutiveMultiFieldStyle()
Y = MultiFieldFESpace(Vector{Float64},[V,Q],mfs)
X = MultiFieldFESpace(Vector{Float64},[U,P],mfs)

@test num_free_dofs(X) == num_free_dofs(U) + num_free_dofs(P)
@test num_free_dofs(X) == num_free_dofs(Y)
@test length(X) == 2
@test num_fields(X) == 2
@test num_fields(V) == 1

dy = get_fe_basis(Y)
dv, dq = dy

dx = get_trial_fe_basis(X)
du, dp = dx

change_domain(dv,PhysicalDomain())

cellmat = integrate(dv*du,quad)
cellvec = integrate(dv*2,quad)
cellmatvec = pair_arrays(cellmat,cellvec)
@test isa(cellmat[end],ArrayBlock)
@test cellmat[1][1,1] != nothing
@test cellmat[1][1,2] == nothing
@test isa(cellvec[end], ArrayBlock)
@test cellvec[1][1] != nothing
@test cellvec[1][2] == nothing

free_values = rand(num_free_dofs(X))
xh = FEFunction(X,free_values)
test_fe_function(xh)
@test isa(xh,FEFunction)
uh, ph = xh
@test isa(uh,FEFunction)
@test isa(ph,FEFunction)

dir_values = zero_dirichlet_values(Y)
@test all(map((dv,Yi) -> dv == zero_dirichlet_values(Yi),dir_values,Y))

xh = FEFunction(X,free_values,dir_values)
@test isa(xh,FEFunction)

cell_isconstr = get_cell_isconstrained(X)
cell_isconstr = get_cell_isconstrained(X,trian)
@test cell_isconstr == Fill(false,num_cells(model))

cell_constr = get_cell_constraints(X)
cell_constr = get_cell_constraints(X,trian)
@test isa(cell_constr,LazyArray{<:Fill{<:BlockMap}})

cell_dof_ids = get_cell_dof_ids(X)
cell_dof_ids = get_cell_dof_ids(X,trian)
@test isa(cell_dof_ids,LazyArray{<:Fill{<:BlockMap}})

cell_is_dirichlet = get_cell_is_dirichlet(X)
cell_is_dirichlet = get_cell_is_dirichlet(X,trian)
@test isa(cell_is_dirichlet,AbstractArray{<:Bool})

cf = CellField(X,get_cell_dof_ids(X,trian))
@test isa(cf,MultiFieldCellField)

test_fe_space(X,cellmatvec,cellmat,cellvec,trian)
test_fe_space(Y,cellmatvec,cellmat,cellvec,trian)

f(x) = sin(4*pi*(x[1]-x[2]^2))+1
fh = interpolate([f,f],X)
fh = interpolate_everywhere([f,f],X)
fh = interpolate_dirichlet([f,f],X)

# BlockMultiFieldStyle

mfs = BlockMultiFieldStyle()
Y = MultiFieldFESpace([V,Q],style=mfs)
X = MultiFieldFESpace([U,P],style=mfs)

fh = interpolate([f,f],X)

x = zero_free_values(X)
interpolate!([f,f],x,X)
@test isa(x,BlockVector)
@test x == get_free_dof_values(fh)

dv = zero_dirichlet_values(X)
interpolate_everywhere!([f,f],x,dv,X)
@test x == get_free_dof_values(fh)

@testset "BlockMultiFieldStyle with non-identity permutation" begin
  model2 = CartesianDiscreteModel((0,1,0,1),(2,2))
  V2 = TestFESpace(model2,ReferenceFE(lagrangian,Float64,2);conformity=:H1)
  Q2 = TestFESpace(model2,ReferenceFE(lagrangian,Float64,1);conformity=:L2)
  U2 = TrialFESpace(V2)
  P2 = TrialFESpace(Q2)
  dΩ2 = Measure(Triangulation(model2),2)

  a2((u,p),(v,q)) = ∫( ∇(v)⋅∇(u) + p*q )dΩ2
  l2((v,q)) = ∫( v + q )dΩ2

  mfs_ref = ConsecutiveMultiFieldStyle()
  X_ref = MultiFieldFESpace([U2,P2];style=mfs_ref)
  Y_ref = MultiFieldFESpace([V2,Q2];style=mfs_ref)
  xh_ref = solve(AffineFEOperator(a2,l2,X_ref,Y_ref))

  mfs_perm = BlockMultiFieldStyle(2,(1,1),(2,1))
  X_perm = MultiFieldFESpace([U2,P2];style=mfs_perm)
  Y_perm = MultiFieldFESpace([V2,Q2];style=mfs_perm)
  xh_perm = solve(AffineFEOperator(a2,l2,X_perm,Y_perm))

  uh_ref, ph_ref = xh_ref
  uh_perm, ph_perm = xh_perm
  @test sum(∫((uh_ref - uh_perm)*(uh_ref - uh_perm))dΩ2) < 1e-24
  @test sum(∫((ph_ref - ph_perm)*(ph_ref - ph_perm))dΩ2) < 1e-24
end

end # module
