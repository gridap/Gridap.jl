
using Test
using Gridap
using Gridap.ReferenceFEs, Gridap.Geometry
using Gridap.Helpers
using Gridap.FESpaces

using LinearAlgebra

p = TRI
D = num_dims(p)
order = 1

reffe = ReferenceFEs.HellanHerrmannJhonsonRefFE(Float64,p,order)
dofs = get_dof_basis(reffe)
prebasis  = get_prebasis(reffe)
shapefuns = get_shapefuns(reffe) 

M = evaluate(dofs, shapefuns)
@test M ≈ I(num_dofs(reffe))

model = simplexify(CartesianDiscreteModel((0,1,0,1), (1,1));positive=false)
Ω = Triangulation(model)
Γ = Boundary(model)
Λ = Skeleton(model)

dΩ = Measure(Ω,4*order)
dΓ = Measure(Γ,4*order)
dΛ = Measure(Λ,4*order)

V = FESpace(model,reffe)
get_cell_dof_ids(V)

#u(x) = SymTensorValue((x[1]*(x[1]-1),0.0,x[2]*(x[2]-1)))
g(x) = VectorValue(x[1]^2+x[2],x[2]^2)
u(x) = TensorValues.symmetric_part((ε(g))(x))
uh = interpolate(u,V)

eh = uh - u
sum(∫(eh⊙eh)dΩ)
collect(get_array(∫(eh⊙eh)dΩ))

nΓ = get_normal_vector(Γ)
collect(get_array(∫(nΓ⋅uh⋅nΓ)dΓ))

ucf = uh# CellField(u,Ω)
nΛ = get_normal_vector(Λ)
get_array(∫((nΛ⋅(ucf⋅nΛ)).plus - (nΛ⋅(ucf⋅nΛ)).minus)dΛ)[1]
get_array(∫((nΛ⋅(ucf⋅nΛ)).plus)dΛ)[1]
get_array(∫((nΛ⋅(ucf⋅nΛ)).minus)dΛ)[1]

pts = get_cell_points(dΩ.quad)
dofs = get_fe_dof_basis(V)
basis = get_fe_basis(V)
dofs(uh)[2]

using Gridap.CellData: SkeletonPair
function LinearAlgebra.dot(a::SkeletonPair,b::SkeletonPair)
  SkeletonPair(dot(a.plus,b.plus),dot(a.minus,b.minus))
end

collect((Operation(.*)(eh,eh))(get_cell_points(dΩ.quad)))[1]

x1 = collect((uh)(get_cell_points(dΩ.quad)))[1]
x2 = collect((ucf)(get_cell_points(dΩ.quad)))[1]
x1 == x2

x3 = collect((eh)(get_cell_points(dΩ.quad)))[1]

inner(x1[1],x1[1])
inner(x1[1],TensorValues.symmetric_part(x1[1]))
