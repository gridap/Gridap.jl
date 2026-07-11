using Gridap, Gridap.MultiField
using Test

# Single field
model = CartesianDiscreteModel((0,1,0,1),(8,8))
order = 1
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

j(u) = ∫(u + im*u*conj(u))dΩ + ∫(im*u*u + 1)dΓ + ∫(mean(u) + im*mean(u*u))dΛ
dj(v,u) = ∫(v + im*v*conj(u) + im*u*conj(v))dΩ + ∫(2im*v*u)dΓ + ∫(mean(v) + 2im*mean(v*u))dΛ

fi = [x-> x[1]*x[2], x->im*x[1]*x[2], x->x[1]+im*x[1]*x[2]]
for f in fi
  for V in [V1,V2]
    xh = interpolate(f,V)
    dj_ad = gradient(j,xh)
    dj_vec = assemble_vector(dj_ad,V)
    dj_vec_analytic = assemble_vector(v -> dj(v,xh),V)
    @test dj_vec≈dj_vec_analytic
  end
end

r(u,v) = ∫(u*v + im*u*conj(u)*v)dΩ + ∫(im*u*u*v + 1*v)dΓ + ∫(mean(u)*mean(v) + im*mean(u*u)*mean(v))dΛ
dr(du,u,v) = ∫(du*v + im*du*conj(u)*v + im*u*conj(du)*v)dΩ + ∫(2im*du*u*v)dΓ + ∫(mean(du)*mean(v) + 2im*mean(du*u)*mean(v))dΛ

for f in fi
  for V in [V1,V2]
    xh = interpolate(f,V)
    dv = get_fe_basis(V)
    dr_ad = jacobian(x->r(x,dv),xh)
    dr_vec = assemble_matrix(dr_ad,V,V)
    dr_vec_analytic = assemble_matrix((u,v) -> dr(u,xh,v),V,V)
    @test dr_vec≈dr_vec_analytic
  end
end

# Multi-field
V1 = FESpace(model,ReferenceFE(lagrangian,Float64,order),vector_type=Vector{ComplexF64})
V2 = FESpace(model, ReferenceFE(lagrangian,Float64,order);vector_type=Vector{ComplexF64},conformity=:L2)
V = MultiFieldFESpace([V1,V2])

j((u1,u2)) = ∫(u1*u2 + im*u1*conj(u1))dΩ + ∫(im*u1*u1*u2 + 1)dΓ + ∫(mean(u1) + im*mean(u1*u1*u2))dΛ
dj((v1,v2),(u1,u2)) = ∫(v1*u2 + v2*u1 + im*v1*conj(u1) + im*u1*conj(v1))dΩ + ∫(2im*v1*u1*u2 + im*u1*u1*v2)dΓ + ∫(mean(v1) + 2im*mean(v1*u1*u2) + im*mean(u1*u1*v2))dΛ

xh = interpolate((x->x[1]*x[2], x->x[1]+im*x[1]*x[2]),V)
dj_ad = gradient(j,xh)
dj_vec = assemble_vector(dj_ad,V)
dj_vec_analytic = assemble_vector(v -> dj(v,xh),V)
@test dj_vec≈dj_vec_analytic

r((u1,u2),(v1,v2)) = ∫(u1*v1 + u2*v2 + im*u1*conj(u1)*v1*u2 + im*u2*conj(u2)*v2*u1)dΩ +
  ∫(im*u1*u1*v1 + im*u2*u2*v2 + 1*v1 + 1*v2)dΓ + ∫(mean(u1)*mean(v1) + mean(u2)*mean(v2) + im*mean(u1*u1)*mean(v1) + im*mean(u2*u2)*mean(v2))dΛ
dr((du1,du2),(u1,u2),(v1,v2)) = ∫(du1*v1 + du2*v2 + im*u1*conj(u1)*v1*du2 + im*u2*conj(u2)*v2*du1 +
  im*du1*conj(u1)*v1*u2 + im*du2*conj(u2)*v2*u1 + im*u1*conj(du1)*v1*u2 + im*u2*conj(du2)*v2*u1)dΩ +
  ∫(2im*du1*u1*v1 + 2im*du2*u2*v2)dΓ + ∫(mean(du1)*mean(v1) + mean(du2)*mean(v2) + 2im*mean(du1*u1)*mean(v1) + 2im*mean(du2*u2)*mean(v2))dΛ

dv = get_fe_basis(V)
dr_ad = jacobian(x->r(x,dv),xh)
dr_vec = assemble_matrix(dr_ad,V,V)
dr_vec_analytic = assemble_matrix((u,v) -> dr(u,xh,v),V,V)
@test dr_vec≈dr_vec_analytic