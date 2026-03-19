module PolytopalFESpacesTests
using Test

using Gridap
using Gridap.Geometry, Gridap.ReferenceFEs, Gridap.TensorValues
using Gridap.Fields, Gridap.FESpaces, Gridap.CellData

using FillArrays

function test_l2_proj(model,V,order,u_exact)
  U = TrialFESpace(V,u_exact)

  Ω = get_triangulation(V)
  dΩ = Measure(Ω,2*order)

  mass_lhs(u,v) = ∫(u⋅v)dΩ
  mass_rhs(v) = ∫(v⋅u_exact)dΩ
  
  op = AffineFEOperator(mass_lhs,mass_rhs,U,V)
  uh = solve(op)
  
  eh = uh - u_exact
  @test sum(∫(eh⋅eh)dΩ) < 1e-10

  uh_i = interpolate(u_exact,U)
  eh = uh - uh_i
  @test sum(∫(eh⋅eh)dΩ) < 1e-10

  d = tempdir()
  writevtk(
    Ω,d*"/polytopal_l2_proj",
    cellfields=["uh" => uh],
    celldata=["cell_id"=>collect(1:num_cells(Ω))],
  )
end

function test_dg_lap(model,V,order,u_exact)
  U = TrialFESpace(V,u_exact)

  Ω = Triangulation(model)
  Γ = Boundary(model)
  Λ = Skeleton(model)
  
  dΩ = Measure(Ω,2*order)
  dΓ = Measure(Γ,2*order)
  dΛ = Measure(Λ,2*order)
  
  β = 100.0
  f(x) = -Δ(u_exact)(x)
  nΛ = get_normal_vector(Λ)
  βΛ = CellField(β ./ get_array(∫(1)dΛ),Λ)
  nΓ = get_normal_vector(Γ)
  βΓ = CellField(β ./ get_array(∫(1)dΓ),Γ)
  
  lap_lhs(u,v) = ∫(∇(u)⋅∇(v))dΩ + 
                 ∫(βΛ*jump(u*nΛ)⋅jump(v*nΛ) - mean(∇(u))⋅jump(v*nΛ) - mean(∇(v))⋅jump(u*nΛ))dΛ + 
                 ∫(βΓ*u*v - (∇(u)⋅nΓ)*v - (∇(v)⋅nΓ)*u)dΓ
  lap_rhs(v) = ∫(v*f)dΩ + ∫(βΓ*u_exact*v - (∇(v)⋅nΓ)*u_exact)dΓ
  
  op = AffineFEOperator(lap_lhs,lap_rhs,U,V)
  uh = solve(op)
  
  eh = uh - u_exact
  @test sum(∫(eh⋅eh)dΩ) < 1e-10
end

order = 2

# 2D bulk
u_exact_2d(x) = x[1]^order + x[2]^order
model = CartesianDiscreteModel((0,1,0,1),(2,2))
pmodel = Geometry.voronoi(Geometry.simplexify(model))

V = FESpaces.PolytopalFESpace(pmodel,Float64,order,space=:P)
test_l2_proj(pmodel,V,order,u_exact_2d)
test_dg_lap(pmodel,V,order,u_exact_2d)

V = FESpaces.PolytopalFESpace(pmodel,Float64,order,space=:P,hierarchical=true,orthonormal=true)
test_l2_proj(pmodel,V,order,u_exact_2d)
test_dg_lap(pmodel,V,order,u_exact_2d)

V = FESpaces.PolytopalFESpace(pmodel,Float64,order,space=:P,local_kernel=:constants)

# 2D skeleton
Γ = Triangulation(ReferenceFE{1},model)

VΓ = FESpaces.PolytopalFESpace(Γ,Float64,order,space=:P,dirichlet_tags=["boundary"])
@test any(get_cell_dof_ids(VΓ).data .< 0)
test_l2_proj(pmodel,VΓ,order,u_exact_2d)

u_exact_2d_vec(x) = VectorValue(x[1]^order,x[2]^order)
VΓ = FESpaces.PolytopalFESpace(Γ,VectorValue{2,Float64},order,space=:P,dirichlet_tags=["boundary"],hierarchical=true,orthonormal=true)
@test any(get_cell_dof_ids(VΓ).data .< 0)
test_l2_proj(pmodel,VΓ,order,u_exact_2d_vec)

# TODO: Does this make sense? 
# VΓ = FESpaces.PolytopalFESpace(Γ,VectorValue{2,Float64},order,space=:P,dirichlet_tags=["boundary"],dirichlet_masks=[true,false])
# @test any(get_cell_dof_ids(VΓ).data .< 0)
# test_l2_proj(pmodel,VΓ,order,u_exact_2d_vec)

# 3D bulk
u_exact_3d(x) = x[1]^order + x[2]^order + x[3]^order
model = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
pmodel = Geometry.PolytopalDiscreteModel(model)

V = FESpaces.PolytopalFESpace(pmodel,Float64,order,space=:P)
test_l2_proj(pmodel,V,order,u_exact_3d)
test_dg_lap(pmodel,V,order,u_exact_3d)

# 3D skeleton
Γ = Triangulation(ReferenceFE{2},model)

VΓ = FESpaces.PolytopalFESpace(Γ,Float64,order,space=:P,dirichlet_tags=["boundary"])
@test any(get_cell_dof_ids(VΓ).data .< 0)
test_l2_proj(pmodel,VΓ,order,u_exact_3d)

u_exact_3d_vec(x) = VectorValue(x[1]^order,x[2]^order,x[3]^order)
VΓ = FESpaces.PolytopalFESpace(Γ,VectorValue{3,Float64},order,space=:P,dirichlet_tags=["boundary"],hierarchical=true,orthonormal=true)
@test any(get_cell_dof_ids(VΓ).data .< 0)
test_l2_proj(pmodel,VΓ,order,u_exact_3d_vec)

# TODO: Does this make sense?
# VΓ = FESpaces.PolytopalFESpace(Γ,VectorValue{3,Float64},order,space=:P,dirichlet_tags=["boundary"],dirichlet_masks=[true,false,true])
# @test any(get_cell_dof_ids(VΓ).data .< 0)
# test_l2_proj(pmodel,VΓ,order,u_exact_3d_vec)

end