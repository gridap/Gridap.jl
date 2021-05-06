module DomainContributionsTests

using Test
using FillArrays
using Gridap.Helpers
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Arrays
using Gridap.CellData
using Gridap.TensorValues
using Gridap.Geometry

domain = (0,1,0,1)
cells = (2,2)
model = simplexify(CartesianDiscreteModel(domain,cells))

Ω = Triangulation(model)
Γ = BoundaryTriangulation(model)
Λ = SkeletonTriangulation(model)
n_Λ = get_normal_vector(Λ)

degree = 2
dΩ = Measure(Ω,degree)
dΓ = Measure(Γ,degree)
dΛ = Measure(Λ,degree)

v = GenericCellField(get_cell_shapefuns(Ω),Ω,ReferenceDomain())
u = GenericCellField(lazy_map(transpose,get_data(v)),v.trian,v.domain_style)

a = ∫(u*v)*dΩ + ∫(u*v)*dΓ + ∫(∇(u)⋅∇(v))*dΩ
@test num_domains(a) == 2
@test Ω in get_domains(a)
@test Γ in get_domains(a)
@test isa(get_contribution(a,Ω),LazyArray)

a = ∫(1)*dΩ + ∫(1)*dΓ
@test sum(a) ≈ 5
@test sum(2*a) ≈ 10
@test sum(a*2) ≈ 10

u = CellField(x->2*x[1],Ω)
v = CellField(x->3*x[2],Ω)

a = ∫(jump(u))*dΛ
@test sum(a) + 1 ≈ 1

a = ∫( (n_Λ.⁺⋅∇(v.⁻))*jump(n_Λ⋅∇(u)) )*dΛ
@test sum(a) + 1 ≈ 1

quad = Quadrature(duffy,2)
dΩ = Measure(Ω,quad)
s = ∫(1)dΩ
@test sum(s) ≈ 1

end # module
