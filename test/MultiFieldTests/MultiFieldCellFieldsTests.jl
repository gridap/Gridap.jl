module MultiFieldCellFieldsTests

using BlockArrays
using FillArrays
using Gridap.Arrays
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.ReferenceFEs
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
x = get_cell_points(quad)

trian_Γ = SkeletonTriangulation(model)
quad_Γ = CellQuadrature(trian_Γ,degree)
x_Γ = get_cell_points(quad_Γ)

V = TestFESpace(model,ReferenceFE(:Lagrangian,VectorValue{2,Float64},order);conformity=:H1)
Q = TestFESpace(model,ReferenceFE(:Lagrangian,Float64,order-1),conformity=:L2)

U = TrialFESpace(V)
P = TrialFESpace(Q)

Y = MultiFieldFESpace([V,Q])
X = MultiFieldFESpace([U,P])

dv, dq = get_cell_shapefuns(Y)
du, dp = get_cell_shapefuns_trial(X)

n = VectorValue(1,2)

cellmat = integrate( (n⋅dv)*dp + dq*dp, quad)
cellvec = integrate( n⋅dv, quad)
@test isa(cellmat,LazyArray{<:Fill{BlockArrayCooMap{2}}})
@test isa(cellvec,LazyArray{<:Fill{BlockArrayCooMap{1}}})

cellmat1 = integrate( ((n⋅dv) - dq)*((n⋅du) + dp), quad)
cellmat2 = integrate( (n⋅dv)*(n⋅du) + (n⋅dv)*dp - dq*(n⋅du) - dq*dp, quad)
test_array(cellmat1,cellmat2,≈)

cellmat1 = integrate( (n⋅dv)*2, quad)
cellmat2 = integrate( (n⋅dv)*fill(2,num_cells(trian)), quad)
test_array(cellmat1,cellmat2,≈)

α = CellField(2,trian)
op(u,∇u,v,∇v,α) = α*(u⋅v) + ∇u⊙∇v
cellmat1 = integrate( op∘(du,∇(du),dv,∇(dv),α) , quad)
cellmat2 = integrate( α*(du⋅dv) + ∇(du)⊙∇(dv) , quad)
test_array(cellmat1,cellmat2,≈)

@test isa(cellmat2,LazyArray{<:Fill{BlockArrayCooMap{2}}})
#@test_broken isa(cellmat1,LazyArray{<:Fill{BlockArrayCooMap{2}}})

α = CellField(2,trian)
op2(u,∇u,α) = α*(∇u⋅u)
cellmat1 = integrate( dv⋅(op2∘(du,∇(du),α)),quad)
cellmat2 = integrate( dv⋅(α*(∇(du)⋅du)),quad)
test_array(cellmat1,cellmat2,≈)
@test isa(cellmat2,LazyArray{<:Fill{BlockArrayCooMap{2}}})
#@test_broken isa(cellmat1,LazyArray{<:Fill{BlockArrayCooMap{2}}})

conv(u,∇u,α) = α*(u⋅∇u)
dconv(du,∇du,u,∇u,α) = conv(u,∇du,α)+conv(du,∇u,α)

u = zero(U)
cellvec2 = integrate(dv⊙(α*(u⋅∇(u))),quad)
cellvec1 = integrate(dv⊙(conv∘(u,∇(u),α)),quad)
test_array(cellvec1,cellvec2,≈)
@test isa(cellvec2,LazyArray{<:Fill{BlockArrayCooMap{1}}})
@test isa(cellvec1,LazyArray{<:Fill{BlockArrayCooMap{1}}})

cellmat1 = integrate( dv⋅(dconv∘(du,∇(du),u,∇(u),α)) , quad)
cellmat2 = integrate( dv⋅( α*(du⋅∇(u)) + α*(u⋅∇(du))), quad)
test_array(cellmat1,cellmat2,≈)
@test isa(cellmat2,LazyArray{<:Fill{BlockArrayCooMap{2}}})
#@test_broken isa(cellmat1,LazyArray{<:Fill{BlockArrayCooMap{2}}})

cellmat_Γ = integrate(  jump(n⋅dv)*dp.⁺ + mean(dq)*jump(dp), quad_Γ)
cellvec_Γ = integrate(  jump(n⋅dv) + mean(dq), quad_Γ)
L = 1
R = 2
@test isa(cellmat_Γ,LazyArray{<:Fill{BlockArrayCooMap{2}}})
@test isa(cellvec_Γ,LazyArray{<:Fill{BlockArrayCooMap{1}}})

cell = 1
@test isa(cellmat_Γ[cell][Block(L,R)],BlockArrayCoo)
@test isa(cellvec_Γ[cell][Block(L)],BlockArrayCoo)

cellmat1_Γ = integrate(((n⋅dv.⁺)-dq.⁻)*((n⋅du.⁺)+dp.⁻),quad_Γ)
cellmat2_Γ = integrate((n⋅dv.⁺)*(n⋅du.⁺)+(n⋅dv.⁺)*dp.⁻-dq.⁻*(n⋅du.⁺)-dq.⁻*dp.⁻,quad_Γ)
test_array(cellmat1_Γ,cellmat2_Γ,≈)


#a = cellmat_Γ
#using BenchmarkTools
#cache = array_cache(a)
#@btime getindex!($cache,$a,2)


end # module
