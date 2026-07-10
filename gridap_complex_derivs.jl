using Gridap, Gridap.MultiField
using Test


model = CartesianDiscreteModel((0,1,0,1),(8,8))
order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
reffe = ReferenceFE(lagrangian,Float64,order)
Ω = Triangulation(model)
dΩ = Measure(Ω,2order)
Γ = Boundary(model)
dΓ = Measure(Γ,2order)
Λ = SkeletonTriangulation(model)
dΛ = Measure(Λ,2order)

V1 = FESpace(model,ReferenceFE(lagrangian,Float64,order);
  vector_type=Vector{ComplexF64})
V2 = FESpace(model, ReferenceFE(lagrangian,Float64,order);
  vector_type=Vector{ComplexF64},conformity=:L2) 

xh1 = interpolate(x -> x[1]*x[2],V1)
xh2 = interpolate(x -> im*x[1]*x[2],V1)
xh3 = interpolate(x -> x[1]+im*x[1]*x[2],V1)
xh3_L2 = interpolate(x -> x[1]+im*x[1]*x[2],V2)

j(u) = ∫(u + im*u*conj(u))dΩ + ∫(im*u*u + 1)dΓ + ∫(mean(u) + im*mean(u*u))dΛ
dj(v,u) = ∫(v + im*v*conj(u) + im*u*conj(v))dΩ + ∫(2im*v*u)dΓ + ∫(mean(v) + 2im*mean(v*u))dΛ

_dj = gradient(j,xh1)
_dj_vec = assemble_vector(_dj,V1)
_dj_vec_analytic = assemble_vector(v -> dj(v,xh1),V1)
_dj_vec≈_dj_vec_analytic

##

_dj = gradient(u->∫(u)dΩ,xh1)
_dj_vec = assemble_vector(_dj,V1)
_dj_vec_analytic = assemble_vector(v -> ∫(v)dΩ,V1)
_dj_vec≈_dj_vec_analytic

_dj = gradient(u->∫(u)dΩ,xh1)
_dj_vec = assemble_vector(_dj,V1)
_dj_vec_analytic = assemble_vector(v-> ∫(v)dΩ,V1)
_dj_vec≈_dj_vec_analytic

_dj = gradient(u->∫(u*u)dΩ,xh1)
_dj_vec = assemble_vector(_dj,V1)
_dj_vec_analytic = assemble_vector(v-> ∫(2*xh1*v)dΩ,V1)
_dj_vec≈_dj_vec_analytic

_dj = gradient(u->∫(im*u*u)dΩ,xh1)
_dj_vec = assemble_vector(_dj,V1)
_dj_vec_analytic = assemble_vector(v-> ∫(2im*xh1*v)dΩ,V1)
_dj_vec≈_dj_vec_analytic

_dj = gradient(u->∫(u + im*u*u)dΩ,xh1)
_dj_vec = assemble_vector(_dj,V1)
_dj_vec_analytic = assemble_vector(v-> ∫(v + 2im*xh1*v)dΩ,V1)
_dj_vec≈_dj_vec_analytic

_dj = gradient(u->∫(u)dΩ,xh2)
_dj_vec = assemble_vector(_dj,V1)
_dj_vec_analytic = assemble_vector(v-> ∫(v)dΩ,V1)
_dj_vec≈_dj_vec_analytic

_dj = gradient(u->∫(im*u)dΩ,xh2)
_dj_vec = assemble_vector(_dj,V1)
_dj_vec_analytic = assemble_vector(v-> ∫(im*v)dΩ,V1)
_dj_vec≈_dj_vec_analytic

_dj = gradient(u->∫(u + im*u*u)dΩ,xh3)
_dj_vec = assemble_vector(_dj,V1)
_dj_vec_analytic = assemble_vector(v-> ∫(v + 2im*xh3*v)dΩ,V1)
_dj_vec≈_dj_vec_analytic

_dj = gradient(u->∫(mean(u))dΛ,xh3_L2)
_dj_vec = assemble_vector(_dj,V2)
_dj_vec_analytic = assemble_vector(v -> ∫(mean(v))dΛ,V2)
_dj_vec≈_dj_vec_analytic

dv = get_fe_basis(V)
_dj = jacobian(u->∫(mean(u)*mean(dv))dΛ,xh3_L2)
_dj_vec = assemble_matrix(_dj,V2,V2)
_dj_vec_analytic = assemble_matrix((u,v) -> ∫(mean(u)*mean(v))dΛ,V2,V2)
_dj_vec≈_dj_vec_analytic

# FAILS
_dj = gradient(u->∫(mean(u*u))dΛ,xh3_L2)
_dj_vec = assemble_vector(_dj,V2)
_dj_vec_analytic = assemble_vector(v -> ∫(mean(2*xh3_L2*v))dΛ,V2)
_dj_vec≈_dj_vec_analytic

dv = get_fe_basis(V)
_dj = jacobian(u->∫(mean(u*u)*jump(dv))dΛ,xh3_L2)
_dj_vec = assemble_matrix(_dj,V2,V2)
_dj_vec_analytic = assemble_matrix((u,v) -> ∫(mean(2*xh3_L2*u)*jump(v))dΛ,V2,V2)
_dj_vec≈_dj_vec_analytic