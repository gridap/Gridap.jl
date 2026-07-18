using Gridap, Gridap.MultiField
using LinearAlgebra
using Random
using Test

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
V2 = FESpace(model,ReferenceFE(lagrangian,Float64,order);
  conformity=:L2,vector_type=Vector{ComplexF64})
h = 1.0e-6

fd_gradient(f,u,du) = real((f(u + h*du) - f(u - h*du))/(2h))
fd_residual(r,u,du) = (r(u + h*du) - r(u - h*du))/(2h)
functional_value(f,V,u) = sum(f(FEFunction(V,u)))
residual_vector(r,V,dv,u) = assemble_vector(r(FEFunction(V,u),dv),V)

@testset "Complex gradients" begin
  # gradient is intended for real-valued objectives. This holomorphic example
  # checks the current extension, which returns the Riesz gradient of real(F).
  F_holo(u) = ∫(u*u)dΩ + ∫(u*u)dΓ + ∫(mean(u)*mean(u))dΛ
  dreal_F_holo(q,u) = ∫(2*q*conj(u))dΩ + ∫(2*q*conj(u))dΓ + ∫(2*mean(q)*conj(mean(u)))dΛ
  J_nonholo(u) = ∫(conj(u)*u)dΩ + ∫(conj(u)*u)dΓ + ∫(mean(conj(u))*mean(u))dΛ
  dJ_nonholo(q,u) = ∫(2*q'*u)dΩ + ∫(2*q'*u)dΓ + ∫(2*mean(q)'*mean(u))dΛ

  for (index,V) in enumerate((V1,V2))
    uh = interpolate(x -> x[1] + im*x[2],V)
    u = get_free_dof_values(uh)
    du = randn(ComplexF64,length(u))

    ad_holo = assemble_vector(gradient(F_holo,uh),V)
    analytic_holo = assemble_vector(q -> dreal_F_holo(q,uh),V)
    @test ad_holo ≈ analytic_holo
    @test real(dot(ad_holo,du)) ≈ fd_gradient(u -> real(functional_value(F_holo,V,u)),u,du) rtol=1.0e-6

    ad_nonholo = assemble_vector(gradient(J_nonholo,uh),V)
    analytic_nonholo = assemble_vector(q -> dJ_nonholo(q,uh),V)
    @test ad_nonholo ≈ analytic_nonholo
    @test real(dot(ad_nonholo,du)) ≈ fd_gradient(u -> real(functional_value(J_nonholo,V,u)),u,du) rtol=1.0e-6
  end
end

@testset "Complex Jacobians" begin
  R_holo(u,v) = ∫(u*u*v)dΩ + ∫(u*u*v)dΓ + ∫(mean(u)*mean(u)*mean(v))dΛ
  dR_holo(du,u,v) = ∫(2*u*du*v)dΩ + ∫(2*u*du*v)dΓ + ∫(2*mean(u)*mean(du)*mean(v))dΛ
  R_nonholo(u,v) = ∫(conj(u)*v)dΩ + ∫(conj(u)*v)dΓ + ∫(mean(conj(u))*mean(v))dΛ
  dR_nonholo(du,u,v) = ∫(conj(du)*v)dΩ + ∫(conj(du)*v)dΓ + ∫(mean(conj(du))*mean(v))dΛ

  for (index,V) in enumerate((V1,V2))
    uh = interpolate(x -> x[1] + im*x[2],V)
    u = get_free_dof_values(uh)
    du = randn(ComplexF64,length(u))
    dv = get_fe_basis(V)

    ad_holo = assemble_matrix(jacobian(u -> R_holo(u,dv),uh),V,V)
    analytic_holo = assemble_matrix((u,v) -> dR_holo(u,uh,v),V,V)
    @test ad_holo ≈ analytic_holo
    @test ad_holo*du ≈ fd_residual(u -> residual_vector(R_holo,V,dv,u),u,du) rtol=1.0e-6

    ad_nonholo = assemble_matrix(jacobian(u -> R_nonholo(u,dv),uh),V,V)
    analytic_real_direction = assemble_matrix((u,v) -> dR_nonholo(u,uh,v),V,V)
    analytic_nonholo_action = assemble_vector(v -> dR_nonholo(FEFunction(V,du),uh,v),V)
    fd_nonholo = fd_residual(u -> residual_vector(R_nonholo,V,dv,u),u,du)
    @test ad_nonholo ≈ analytic_real_direction
    @test analytic_nonholo_action ≈ fd_nonholo rtol=1.0e-6
    @test_broken ad_nonholo*du ≈ fd_nonholo rtol=1.0e-6
  end
end

@testset "Complex multifield AD" begin
  for (style_index,style) in enumerate((ConsecutiveMultiFieldStyle(),BlockMultiFieldStyle()))
    V = MultiFieldFESpace([V1,V2];style)
    uh = interpolate((x -> x[1] + im*x[2],x -> x[1] - im*x[2]),V)
    u = get_free_dof_values(uh)
    du = randn(ComplexF64,length(u))
    dv = get_fe_basis(V)

    F_holo((u1,u2)) = ∫(u1*u2)dΩ + ∫(u1*u2)dΓ + ∫(mean(u1)*mean(u2))dΛ
    dreal_F_holo((q1,q2),(u1,u2)) = ∫(q1*conj(u2) + q2*conj(u1))dΩ + ∫(q1*conj(u2) + q2*conj(u1))dΓ + ∫(mean(q1)*conj(mean(u2)) + mean(q2)*conj(mean(u1)))dΛ
    J_nonholo((u1,u2)) = ∫(conj(u1)*u1 + conj(u2)*u2)dΩ + ∫(conj(u1)*u1 + conj(u2)*u2)dΓ + ∫(mean(conj(u1))*mean(u1) + mean(conj(u2))*mean(u2))dΛ
    dJ_nonholo((q1,q2),(u1,u2)) = ∫(2*q1'*u1 + 2*q2'*u2)dΩ + ∫(2*q1'*u1 + 2*q2'*u2)dΓ + ∫(2*mean(q1)'*mean(u1) + 2*mean(q2)'*mean(u2))dΛ

    for ad_type in (:monolithic,:split)
      ad_holo = assemble_vector(gradient(F_holo,uh;ad_type),V)
      analytic_holo = assemble_vector(q -> dreal_F_holo(q,uh),V)
      @test ad_holo ≈ analytic_holo
      @test real(dot(ad_holo,du)) ≈ fd_gradient(u -> real(functional_value(F_holo,V,u)),u,du) rtol=1.0e-6

      ad_nonholo = assemble_vector(gradient(J_nonholo,uh;ad_type),V)
      analytic_nonholo = assemble_vector(q -> dJ_nonholo(q,uh),V)
      @test ad_nonholo ≈ analytic_nonholo
      @test real(dot(ad_nonholo,du)) ≈ fd_gradient(u -> real(functional_value(J_nonholo,V,u)),u,du) rtol=1.0e-6
    end

    R_holo((u1,u2),(v1,v2)) = ∫(u1*u2*v1 + u2*u2*v2)dΩ + ∫(u1*u2*v1 + u2*u2*v2)dΓ + ∫(mean(u1)*mean(u2)*mean(v1) + mean(u2)*mean(u2)*mean(v2))dΛ
    dR_holo((du1,du2),(u1,u2),(v1,v2)) = ∫((du1*u2 + u1*du2)*v1 + 2*u2*du2*v2)dΩ + ∫((du1*u2 + u1*du2)*v1 + 2*u2*du2*v2)dΓ + ∫((mean(du1)*mean(u2) + mean(u1)*mean(du2))*mean(v1) + 2*mean(u2)*mean(du2)*mean(v2))dΛ
    R_nonholo((u1,u2),(v1,v2)) = ∫(conj(u1)*v1 + conj(u2)*v2)dΩ + ∫(conj(u1)*v1 + conj(u2)*v2)dΓ + ∫(mean(conj(u1))*mean(v1) + mean(conj(u2))*mean(v2))dΛ
    dR_nonholo((du1,du2),(u1,u2),(v1,v2)) = ∫(conj(du1)*v1 + conj(du2)*v2)dΩ + ∫(conj(du1)*v1 + conj(du2)*v2)dΓ + ∫(mean(conj(du1))*mean(v1) + mean(conj(du2))*mean(v2))dΛ

    for ad_type in (:monolithic,:split)
      ad_holo = assemble_matrix(jacobian(u -> R_holo(u,dv),uh;ad_type),V,V)
      analytic_holo = assemble_matrix((u,v) -> dR_holo(u,uh,v),V,V)
      @test ad_holo ≈ analytic_holo
      @test ad_holo*du ≈ fd_residual(u -> residual_vector(R_holo,V,dv,u),u,du) rtol=1.0e-6

      ad_nonholo = assemble_matrix(jacobian(u -> R_nonholo(u,dv),uh;ad_type),V,V)
      analytic_real_direction = assemble_matrix((u,v) -> dR_nonholo(u,uh,v),V,V)
      analytic_nonholo_action = assemble_vector(v -> dR_nonholo(FEFunction(V,du),uh,v),V)
      fd_nonholo = fd_residual(u -> residual_vector(R_nonholo,V,dv,u),u,du)
      @test ad_nonholo ≈ analytic_real_direction
      @test analytic_nonholo_action ≈ fd_nonholo rtol=1.0e-6
      @test_broken ad_nonholo*du ≈ fd_nonholo rtol=1.0e-6
    end
  end
end