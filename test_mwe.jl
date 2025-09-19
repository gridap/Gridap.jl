# using Pkg; Pkg.activate("./Gridap");
using Gridap
using Gridap.Algebra
using Gridap.Arrays
using Gridap.Geometry
using Gridap.CellData
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.MultiField

using Gridap.Arrays: lazy_collect
using FillArrays

using Test

## SingleField
model = CartesianDiscreteModel((0,1,0,1),(3,3))
Ω1 = Triangulation(model,[1,2,3])
Λ = SkeletonTriangulation(model)

V = FESpace(Ω1,ReferenceFE(lagrangian,Float64,1))
uh = zero(V);
fill!(get_free_dof_values(uh),1.0)
dΛ = Measure(Λ,2)

# Gradient
f2(xh) = ∫(mean(xh)*mean(xh))dΛ
dv = get_fe_basis(V);
j = gradient(f2,uh)
J = assemble_vector(j,V)

df2(dxh,xh) = ∫(2*mean(dxh)*mean(xh))dΛ
J_analytic = assemble_vector(dv->df2(dv,uh),V)

@test J ≈ J_analytic

# Jacobian
f2(xh,yh) = ∫(mean(xh)*mean(xh)*mean(yh))dΛ
dv = get_fe_basis(V);
j = jacobian(uh->f2(uh,dv),uh)
J = assemble_matrix(j,V,V)

df2(xh,dxh,yh) = ∫(2*mean(dxh)*mean(xh)*mean(yh))dΛ
op = FEOperator(f2,df2,V,V)
J_analytic = jacobian(op,uh)

@test J ≈ J_analytic

## MultiField
for mod_type in (:skel_trian, :trian, :boundary)
  for style in (ConsecutiveMultiFieldStyle(),BlockMultiFieldStyle())
    model = CartesianDiscreteModel((0,1,0,1),(3,3))
    Ω1 = Triangulation(model,[1,2,3])
    if mod_type == :skel_trian
      Λ = SkeletonTriangulation(model)
    elseif mod_type == :trian
      Λ = Triangulation(model)
    elseif mod_type == :boundary
      Λ = BoundaryTriangulation(model)
    else
      error()
    end

    V1 = FESpace(Ω1,ReferenceFE(lagrangian,Float64,1),conformity=:L2)
    V2 = FESpace(Ω1,ReferenceFE(lagrangian,VectorValue{2,Float64},1),conformity=:L2)
    V3 = FESpace(Ω1,ReferenceFE(lagrangian,Float64,1),conformity=:L2)
    V = MultiFieldFESpace([V1,V2,V3];style)
    uh = zero(V);
    fill!(get_free_dof_values(uh),1.0)
    dΛ = Measure(Λ,2)

    # Gradient
    f(xh) = ∫(mean(xh[1])+mean(xh[2])⋅mean(xh[2])+mean(xh[1])*mean(xh[3]))dΛ
    dv = get_fe_basis(V);
    j = gradient(f,uh;ad_type=:split)
    J = assemble_vector(j,V)

    j_mono = gradient(f,uh;ad_type=:monolithic)
    J_mono = assemble_vector(j_mono,V)

    df2(v,xh) = ∫(mean(v[1])+2*mean(v[2])⋅mean(xh[2])+mean(v[1])*mean(xh[3])+mean(xh[1])*mean(v[3]))dΛ
    J_analytic = assemble_vector(dv->df2(dv,uh),V)
    @test J ≈ J_analytic
    @test J_mono ≈ J_analytic

    # Jacobian
    f2(xh,yh) = ∫(mean(xh[1])⋅mean(yh[1])+mean(xh[2])⋅mean(yh[2])+mean(xh[1])⋅mean(xh[2])⋅mean(yh[2])+mean(xh[1])*mean(xh[3])*mean(yh[3]))dΛ
    dv = get_fe_basis(V);
    j = jacobian(uh->f2(uh,dv),uh;ad_type=:split)
    J = assemble_matrix(j,V,V)

    j_mono = jacobian(uh->f2(uh,dv),uh;ad_type=:monolithic)
    J_mono = assemble_matrix(j_mono,V,V)

    df2(xh,dxh,yh) = ∫(mean(dxh[1])⋅mean(yh[1])+mean(dxh[2])⋅mean(yh[2])+mean(dxh[1])⋅mean(xh[2])⋅mean(yh[2]) +
      mean(xh[1])⋅mean(dxh[2])⋅mean(yh[2])+mean(dxh[1])*mean(xh[3])*mean(yh[3])+mean(xh[1])*mean(dxh[3])*mean(yh[3]))dΛ
    op = FEOperator(f2,df2,V,V)
    J_analytic = jacobian(op,uh)

    @test J ≈ J_analytic
  end
end

#########

## MultiField (no test without skeleton)
# for style in (ConsecutiveMultiFieldStyle(),BlockMultiFieldStyle())
model = CartesianDiscreteModel((0,1,0,1),(3,3))
Ω1 = Triangulation(model,[1,2,3])
Ω2 = Triangulation(model)

V1 = FESpace(Ω1,ReferenceFE(lagrangian,Float64,1))
V2 = FESpace(Ω1,ReferenceFE(lagrangian,VectorValue{2,Float64},1))
V3 = FESpace(Ω1,ReferenceFE(lagrangian,Float64,1))
V = MultiFieldFESpace([V1,V2,V3])#;style)
uh = zero(V);
fill!(get_free_dof_values(uh),1.0)
dΛ = Measure(Ω2,2)

# Gradient
f2(xh) = ∫((xh[1])+(xh[2])⋅(xh[2])+(xh[1])*(xh[3]))dΛ
dv = get_fe_basis(V);
j = gradient(f2,uh;ad_type=:split)
J = assemble_vector(j,V)

j_mono = gradient(f2,uh;ad_type=:monolithic)
J_mono = assemble_vector(j_mono,V)

df2(v,xh) = ∫((v[1])+2*(v[2])⋅(xh[2])+(v[1])*(xh[3])+(xh[1])*(v[3]))dΛ
J_analytic = assemble_vector(dv->df2(dv,uh),V)
@test J ≈ J_analytic
@test J_mono ≈ J_analytic

# Jacobian
f2(xh,yh) = ∫((xh[1])⋅(yh[1])+(xh[2])⋅(yh[2])+(xh[1])⋅(xh[2])⋅(yh[2])+(xh[1])*(xh[3])*(yh[3]))dΛ
dv = get_fe_basis(V);
j = jacobian(uh->f2(uh,dv),uh;ad_type=:split)
J = assemble_matrix(j,V,V)

j_mono = jacobian(uh->f2(uh,dv),uh;ad_type=:monolithic)
J_mono = assemble_matrix(j_mono,V,V)

df2(xh,dxh,yh) = ∫((dxh[1])⋅(yh[1])+(dxh[2])⋅(yh[2])+(dxh[1])⋅(xh[2])⋅(yh[2]) +
  (xh[1])⋅(dxh[2])⋅(yh[2])+(dxh[1])*(xh[3])*(yh[3])+(xh[1])*(dxh[3])*(yh[3]))dΛ
op = FEOperator(f2,df2,V,V)
J_analytic = jacobian(op,uh)

@test J ≈ J_analytic