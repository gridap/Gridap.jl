
using BenchmarkTools
using SparseArrays, LinearAlgebra
using Gridap, Gridap.FESpaces, Gridap.Algebra
using Gridap.CellData, Gridap.Arrays, Gridap.Fields

model = CartesianDiscreteModel((0,1,0,1),(2,2))

order = 1
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffe)
X = MultiFieldFESpace([V,V])

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
Λ = SkeletonTriangulation(model)

degree = 2*order
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
dΛ = Measure(Λ,degree)

n_Γ = get_normal_vector(Γ)
n_Λ = get_normal_vector(Λ)
γ = 1.0
γ0 = 1.0
h = 1.0
f = 1.0
g = 10.0

a((u,p),(v,q)) =
  ∫( ∇(v)⊙∇(u) - ∇(q)⋅u + v⋅∇(p) )*dΩ +
  ∫( (γ/h)*v⋅u - v⋅(n_Γ⋅∇(u)) - (n_Γ⋅∇(v))⋅u + 2*(q*n_Γ)⋅u )*dΓ +
  ∫(
    (γ/h)*jump(v⊗n_Λ)⊙jump(u⊗n_Λ) -
      jump(v⊗n_Λ)⊙mean(∇(u)) -
      mean(∇(v))⊙jump(u⊗n_Λ)  +
      (γ0*h)*jump(q*n_Λ)⋅jump(p*n_Λ) +
      jump(q*n_Λ)⋅mean(u) -
      mean(v)⋅jump(p*n_Λ)
   )*dΛ

l((v,q)) =
  ∫( v⋅f + q*g )*dΩ +
  ∫( (γ/h)*v⋅u - (n_Γ⋅∇(v))⋅u + (q*n_Γ)⋅u )*dΓ

op = AffineFEOperator(a,l,X,X)


u, p = get_trial_fe_basis(X)
v, q = get_fe_basis(X)

u = get_trial_fe_basis(V)
v = get_fe_basis(V)

cf1 = jump(v⊗n_Λ)
cf2 = jump(u⊗n_Λ)
cf3 = cf1 ⊙ cf2

cf4 = mean(∇(u))
cf5 = cf1 ⊙ cf4

cf1_data = testitem(get_data(cf1))
cf4_data = testitem(get_data(cf4))

_cf4 = change_domain(cf4,Λ,ReferenceDomain())
_cf4_data = testitem(get_data(_cf4))
return_value(Broadcasting(Operation(inner)),cf1_data,_cf4_data)

get_data(cf5)

args = map(get_data,cf5.args)
lazy_map(Broadcasting(cf5.op),args...)

pts = CellData._get_cell_points(cf1,cf2)
x = testitem(get_data(pts))
f1 = testitem(get_data(cf1))
f2 = testitem(get_data(cf2))

fx1 = return_value(f1,x)
fx2 = return_value(f2,x)
r = Fields.BroadcastingFieldOpMap(⊙)(fx1,fx2)

cf4 = v⊗n_Λ
cf4_p = testitem(get_data(cf4.plus))
cf4_m = testitem(get_data(cf4.minus))

get_data(jump(n_Λ))

fkx = map(cfk -> evaluate(testitem(get_data(cfk)),x),cf1.args)
r = Fields.BroadcastingFieldOpMap(+)(fkx...)

