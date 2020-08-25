module MultiFieldCellFieldsTests

using BlockArrays
using Gridap.Arrays
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Fields
using Gridap.Integration
using Gridap.CellData
using Gridap.MultiField
using Gridap.TensorValues
using Test

domain = (0,1,0,1)
cells = (2,2)
model = CartesianDiscreteModel(domain,cells)

trian = Triangulation(model)

u1(x) = sin(x[1])
cf1 = CellField(u1,trian)

u2(x) = cos(x[2])
cf2 = CellField(u2,trian)

cf = MultiFieldCellField([cf1,cf2])

@test cf1 === cf[1]
@test cf2 === cf[2]

_cf1, _cf2 = cf

@test cf1 === _cf1
@test cf2 === _cf2

order = 2

domain = (0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

degree = order

trian = get_triangulation(model)
quad = CellQuadrature(trian,degree)
q = get_coordinates(quad)

trian_Γ = SkeletonTriangulation(model)
quad_Γ = CellQuadrature(trian_Γ,degree)
q_Γ = get_coordinates(quad_Γ)

V = TestFESpace(model=model,order=order,reffe=:Lagrangian,conformity=:H1,valuetype=VectorValue{2,Float64})
Q = TestFESpace(model=model,order=order-1,reffe=:Lagrangian,conformity=:L2,valuetype=Float64)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

dv, dq = get_cell_basis(Y)
du, dp = get_cell_basis(X)

n = VectorValue(1,2)

cellmat = integrate( (n⋅dv)*dp + dq*dp, trian, quad)
cellvec = integrate( n⋅dv, trian, quad)
@test isa(cellmat,VectorOfBlockArrayCoo)
@test isa(cellvec,VectorOfBlockArrayCoo)

cellmat1 = integrate( ((n⋅dv) - dq)*((n⋅du) + dp), trian, quad)
cellmat2 = integrate( (n⋅dv)*(n⋅du) + (n⋅dv)*dp - dq*(n⋅du) - dq*dp, trian, quad)
test_array(cellmat1,cellmat2,≈)

cellmat1 = integrate( (n⋅dv)*2, trian, quad)
cellmat2 = integrate( (n⋅dv)*fill(2,num_cells(trian)), trian, quad)
test_array(cellmat1,cellmat2,≈)

op(u,∇u,v,∇v) = u⋅v + ∇u⊙∇v

cellmat1 = integrate( operate(op,du,∇(du),dv,∇(dv)) , trian, quad)
cellmat2 = integrate( du⋅dv + ∇(du)⊙∇(dv) , trian, quad)
test_array(cellmat1,cellmat2,≈)

dv_Γ, dq_Γ = restrict(get_cell_basis(Y), trian_Γ)
du_Γ, dp_Γ = restrict(get_cell_basis(X), trian_Γ)

cellmat_Γ = integrate(  jump(n⋅dv_Γ)*dp_Γ.⁺ + mean(dq_Γ)*jump(dp_Γ), trian_Γ,quad_Γ)
cellvec_Γ = integrate(  jump(n⋅dv_Γ) + mean(dq_Γ), trian_Γ,quad_Γ)
L = 1
R = 2
@test isa(cellmat_Γ,VectorOfBlockArrayCoo)
@test isa(cellmat_Γ[Block(L,R)],VectorOfBlockArrayCoo)
@test isa(cellvec_Γ,VectorOfBlockArrayCoo)
@test isa(cellvec_Γ[Block(L)],VectorOfBlockArrayCoo)

cell = 1
@test isa(cellmat_Γ[cell][Block(L,R)],BlockArrayCoo)
@test isa(cellvec_Γ[cell][Block(L)],BlockArrayCoo)

cellmat1_Γ = integrate(((n⋅dv_Γ.⁺)-dq_Γ.⁻)*((n⋅du_Γ.⁺)+dp_Γ.⁻),trian_Γ,quad_Γ)
cellmat2_Γ = integrate((n⋅dv_Γ.⁺)*(n⋅du_Γ.⁺)+(n⋅dv_Γ.⁺)*dp_Γ.⁻-dq_Γ.⁻*(n⋅du_Γ.⁺)-dq_Γ.⁻*dp_Γ.⁻,trian_Γ,quad_Γ)
test_array(cellmat1_Γ,cellmat2_Γ,≈)


#a = cellmat_Γ
#using BenchmarkTools
#cache = array_cache(a)
#@btime getindex!($cache,$a,2)


end # module
