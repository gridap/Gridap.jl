module UniformRefinementTests
  
using Test
using Gridap
using Gridap.Adaptivity
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.Adaptivity: uniformly_refine

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

end
