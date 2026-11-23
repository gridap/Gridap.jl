module PolytopalFESpacesTests
using Test

using Gridap
using Gridap.Geometry, Gridap.ReferenceFEs, Gridap.TensorValues
using Gridap.Fields, Gridap.FESpaces, Gridap.CellData
using Gridap.Arrays, Gridap.Polynomials

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

# PatchFESpaces

model = CartesianDiscreteModel((0,1,0,1),(2,2))
ptopo = Geometry.PatchTopology(get_grid_topology(model), Table([[1,2],[3,4]]))

pspace = FESpaces.PatchFESpace(model,ptopo,Float64,order,space=:P)
pspace = FESpaces.PatchFESpace(model,ptopo,Float64,order,space=:P,idegree=2*order+2)
@test pspace.idegree == 2*order+2

# Interpolation degree: Projecting P_2 data onto P_0 requires idegree = 2 to get exact means.

model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(2,2)))
pmodel = Gridap.Geometry.voronoi(Gridap.Geometry.simplexify(model))
Ω = Triangulation(pmodel)
dΩ = Measure(Ω,4)
p2(x) = x[1]^2 + x[1]*x[2]

V0 = FESpaces.PolytopalFESpace(pmodel,Float64,0,space=:P)
@test V0.idegree == 0
V0 = FESpaces.PolytopalFESpace(pmodel,Float64,0,space=:P,idegree=2)
@test V0.idegree == 2
uh = interpolate(p2,V0)
@test maximum(abs,get_array(∫(uh - p2)dΩ)) < 1.e-12

VΓ0 = FESpaces.PolytopalFESpace(Triangulation(ReferenceFE{1},pmodel),Float64,0,space=:P,dirichlet_tags=["boundary"],idegree=2)
uh = interpolate_dirichlet(p2,VΓ0)
dΓ = Measure(BoundaryTriangulation(pmodel),4)
@test maximum(abs,get_array(∫(uh - p2)dΓ)) < 1.e-12

# Orthonormal bases: the mass matrix is the identity, on polytopal grids (physical domain),
# standard grids (reference domain) and patch grids (physical domain on a standard grid).

function test_orthonormal(V,order)
  dΩ = Measure(get_triangulation(V),2*order)
  u, v = get_trial_fe_basis(V), get_fe_basis(V)
  for M in get_array(∫(u⋅v)dΩ)
    @test M ≈ one(M)
  end
end

model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,2,0,1),(4,3)))
pmodel = Gridap.Geometry.voronoi(Gridap.Geometry.simplexify(model))
test_orthonormal(FESpaces.PolytopalFESpace(pmodel,Float64,order,orthonormal=true),order)
test_orthonormal(FESpaces.PolytopalFESpace(Triangulation(ReferenceFE{1},pmodel),Float64,order,orthonormal=true),order)
test_orthonormal(FESpaces.PolytopalFESpace(Triangulation(ReferenceFE{1},pmodel),VectorValue{2,Float64},order,orthonormal=true),order)

# The centroid map is isotropic: same scaling in all directions, also for stretched cells.

for cmap in FESpaces.get_centroid_map(get_grid(pmodel))
  G = cmap.gradient
  @test G ≈ G[1,1]*one(G)
end

model = CartesianDiscreteModel((0,2,0,1,0,1),(2,2,2))
pmodel = Geometry.PolytopalDiscreteModel(model)
Γ = Triangulation(ReferenceFE{2},pmodel)
for (cmap,K) in zip(FESpaces.get_centroid_map(get_grid(Γ)),FESpaces.get_tangent_frames(get_grid(Γ)))
  @test transpose(K)⋅K ≈ one(TensorValue{2,2,Float64})
  @test cmap.gradient ≈ (cmap.gradient⊙K/2)*K
end

# Tangent-valued spaces on embedded faces: the values are vectors of the ambient space,
# tangent to the face.

VΓ = FESpaces.PolytopalFESpace(Γ,VectorValue{2,Float64},order,orthonormal=true)
test_orthonormal(VΓ,order)
normal(K) = VectorValue(K[1,1],K[2,1],K[3,1]) × VectorValue(K[1,2],K[2,2],K[3,2])
nΓ = CellField(map(normal,FESpaces.get_tangent_frames(get_grid(Γ))),Γ)
vΓ = get_fe_basis(VΓ)
pts = get_cell_points(Measure(Γ,2*order))
@test eltype(eltype(vΓ(pts))) == VectorValue{3,Float64}
@test all(x -> maximum(abs,x) < 1.e-12, (vΓ⋅nΓ)(pts))

t_exact(x) = VectorValue(x[2]^order,x[3]^order - x[1],x[1]*x[2])
tΓ = t_exact - (t_exact⋅nΓ)*nΓ
eh = interpolate(tΓ,VΓ) - tΓ
@test sum(∫(eh⋅eh)Measure(Γ,2*order)) < 1e-10

# Spaces from a given prebasis

mono = MonomialBasis(Val(2),VectorValue{2,Float64},order,Polynomials._p_filter)
V1 = FESpaces.PolytopalFESpace(Γ,VectorValue{2,Float64},order)
V2 = FESpaces.PolytopalFESpace(Γ,mono)
@test V2.order == order
@test get_fe_basis(V1)(pts) == get_fe_basis(V2)(pts)

n = length(mono)
C = zeros(n,n-2); for i in 1:n-2; C[i,i] = 1.0; end
V3 = FESpaces.PolytopalFESpace(pmodel,linear_combination(C,mono);trian=Γ,orthonormal=true,dirichlet_tags=["boundary"])
@test num_free_dofs(V3) + num_dirichlet_dofs(V3) == (n-2)*num_cells(Γ)
@test num_dirichlet_dofs(V3) > 0
test_orthonormal(V3,order)

end
