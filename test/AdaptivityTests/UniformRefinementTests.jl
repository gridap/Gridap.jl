module UniformRefinementTests
  
using Test
using Gridap
using Gridap.Adaptivity
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Adaptivity: uniformly_refine
using Gridap.CellData: change_domain

using ..EdgeBasedRefinementTests: test_grid_transfers

# Setup base models
has_affine_map = true
cart_model = CartesianDiscreteModel((0,1,0,1),(4,4))
model1 = UnstructuredDiscreteModel(cart_model;has_affine_map)
model2 = simplexify(model1)

cart_model = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
model3 = UnstructuredDiscreteModel(cart_model;has_affine_map)
model4 = simplexify(model3)

model5 = Geometry.DiscreteModelMock()

periodic_model = CartesianDiscreteModel((0,1,0,1),(4,4);isperiodic=(true,false))
model6 = UnstructuredDiscreteModel(periodic_model;has_affine_map)

n = 3

visualize = false
if visualize
  path = mkpath("tmp/")
end

# QUAD
ref_model = refine(model1,n)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_quad_$n"))
test_grid_transfers(model1,ref_model,1)

# TRI
ref_model = refine(model2,n)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_tri_$n"))
test_grid_transfers(model2,ref_model,1)

# HEX
ref_model = refine(model3,n)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_hex_$n"))
test_grid_transfers(model3,ref_model,1)

# TET
ref_model = refine(model4,n)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_tet_$n"))
test_grid_transfers(model4,ref_model,1)

# Mock
ref_model = refine(model5,n)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_mock_$n"))
test_grid_transfers(model5,ref_model,1)

# Periodic QUAD
ref_model = refine(model6,n)
visualize && writevtk(ref_model,joinpath(path,"uniform_periodic_quad_$n"))
# test_grid_transfers(model6,ref_model,1)

# Partial refinement
cell_refine_masks = falses(num_cells(model5))
cell_refine_masks[[1,3,5]] .= true
ref_model = uniformly_refine(model5,n,cell_refine_masks)
visualize && writevtk(Triangulation(ref_model.model),joinpath(path,"uniform_mock_partial_$n"))
test_grid_transfers(model5,ref_model,1)

# Simplex refine(model,2) supports RT and Nedelec
ref_model = refine(model2,2)
Ω = Triangulation(ref_model.model)
reffe = ReferenceFE(TRI,raviart_thomas,1)
v2(x) = VectorValue(-0.5*x[1]+1.0,-0.5*x[2])
V = FESpace(ref_model.model,reffe,conformity=:HDiv)
vh = interpolate(v2,V)
dΩ = Measure(Ω,2)
e = v2 - vh
@test sqrt(sum(∫( e⋅e )*dΩ )) < 1.0e-10
reffe = ReferenceFE(TRI,nedelec,0)
u2((x,y)) = 2*VectorValue(-y,x)
V = TestFESpace(ref_model.model,reffe,dirichlet_tags = "boundary")
U = TrialFESpace(V,u2)
uh = interpolate(u2,U)
dΩ = Measure(Ω,1)
e = u2 - uh
@test sqrt(sum(∫( e⋅e )*dΩ )) < 1.0e-10
ref_model = refine(model4,2)
Ω = Triangulation(ref_model.model)
reffe = ReferenceFE(TET,raviart_thomas,1)
v3(x) = VectorValue(-0.5*x[1]+1.0,-0.5*x[2],-0.5*x[3])
V = FESpace(ref_model.model,reffe,conformity=:HDiv)
vh = interpolate(v3,V)
dΩ = Measure(Ω,2)
e = v3 - vh
@test sqrt(sum(∫( e⋅e )*dΩ )) < 1.0e-10
reffe = ReferenceFE(TET,nedelec,0)
u3((x,y,z)) = 2*VectorValue(-y,x,0.0) - VectorValue(0.0,-z,y)
V = TestFESpace(ref_model.model,reffe,dirichlet_tags = "boundary")
U = TrialFESpace(V,u3)
uh = interpolate(u3,U)
dΩ = Measure(Ω,1)
e = u3 - uh
@test sqrt(sum(∫( e⋅e )*dΩ )) < 1.0e-10

# Surface mesh in 3D: the boundary of a tetrahedron
coords = Point{3,Float64}[(0,0,0),(2,0,0),(0,2,0),(0,0,2)]
cell_vertices = Table([[1,2,3],[1,2,4],[1,3,4],[2,3,4]])
grid = UnstructuredGrid(coords,cell_vertices,[LagrangianRefFE(Float64,TRI,1)],ones(Int8,4))
surface = UnstructuredDiscreteModel(grid)
fine_surface = refine(surface,2).model
midpoints = Point{3,Float64}[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(1,0,1),(0,1,1)]
fine_coords = get_vertex_coordinates(get_grid_topology(fine_surface))
@test fine_coords isa Vector{Point{3,Float64}}
@test issetequal(fine_coords,[coords; midpoints])
@test num_cells(fine_surface) == 16

# The type of the new coordinates is inferred
topo = get_grid_topology(surface)
cell_rr, perms = Adaptivity._uniform_cell_rr_and_face_perms(Val(2),topo,collect(get_cell_type(surface)),trues(4),2)
cell_l2g, n_nodes = Adaptivity._uniform_cell_l2gmap_and_nnodes(topo,cell_rr,perms)
cell_ref_coords = lazy_map(get_node_coordinates,lazy_map(Adaptivity.get_ref_grid,cell_rr))
cell_map = get_cell_map(grid)
@test fine_coords == @inferred Adaptivity._uniform_coordinates(cell_ref_coords,cell_l2g,cell_map,n_nodes)

# CellFields move between the levels of a surface, coarse to fine and back
f(x) = VectorValue(x[1]^2, x[2]*x[3], 1.0)
adapted_surface = refine(surface,2)
Ωc, Ωf = Triangulation(surface), Triangulation(adapted_surface)
dΩc, dΩf = Measure(Ωc,4), Measure(Ωf,4)
for fc in (CellField(f,Ωc), change_domain(CellField(f,Ωc),ReferenceDomain()))
  ef = change_domain(fc,Ωf,ReferenceDomain()) - CellField(f,Ωf)
  @test sqrt(sum(∫( ef⋅ef )*dΩf)) < 1.0e-12
end
ec = change_domain(CellField(f,Ωf),Ωc,ReferenceDomain()) - CellField(f,Ωc)
@test sqrt(sum(∫( ec⋅ec )*dΩc)) < 1.0e-12

end
