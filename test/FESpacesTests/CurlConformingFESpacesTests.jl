module CurlConformingFESpacesTests

using Test
using LinearAlgebra
using FillArrays
using Gridap
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.TensorValues
using Gridap.CellData
using Gridap.Fields
using Gridap.ReferenceFEs

domain =(0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
order = 1
u((x,y)) = 2*VectorValue(2,3)
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
Ω = Triangulation(model)
dΩ = Measure(Ω,order)
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10
#using Gridap.Visualization
#writevtk(Ω,"nedel",nsubcells=10,cellfields=["err"=>e,"u"=>u,"uh"=>uh])

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition) |> simplexify
order = 0
u((x,y)) = 2*VectorValue(-y,x)
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
Ω = Triangulation(model)
dΩ = Measure(Ω,order)
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

order = 1
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
Ω = Triangulation(model)
dΩ = Measure(Ω,order)
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

order = 2
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
Ω = Triangulation(model)
dΩ = Measure(Ω,order)
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition) |> simplexify
Ω = Triangulation(model)
dΩ = Measure(Ω,order)
u((x,y,z)) = 2*VectorValue(-y,x,0.) - VectorValue(0.,-z,y)
order = 0
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10
order = 1
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10
order = 2
reffe = ReferenceFE(nedelec,order)
V = TestFESpace(model,reffe,dirichlet_tags = "boundary")
test_single_field_fe_space(V)
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-10

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)

order = 2

u(x) = VectorValue(x[1]*x[1],x[1]*x[1]*x[1],0.0)
#u(x) = VectorValue(2,3,5)
# u(x) = x

reffe = ReferenceFE(nedelec,order)


V = TestFESpace(model,reffe,dirichlet_tags = [21,22])
# dirichlet_tags = "boundary")
test_single_field_fe_space(V)

U = TrialFESpace(V,u)

uh = interpolate(u,U)

e = u - uh

Ω = Triangulation(model)
dΩ = Measure(Ω,order)

el2 = sqrt(sum( ∫( e⋅e )*dΩ ))
@test el2 < 1.0e-9

domain = (0,1,0,1,0,1)
cells  = (2,2,2)
model  = CartesianDiscreteModel(domain,cells)

# Restrict model to cube surface
labels = get_face_labeling(model)
bgface_to_mask = get_face_mask(labels,"boundary",2)
Γface_to_bgface = findall(bgface_to_mask)
Dc2Dp3model = DiscreteModelPortion(DiscreteModel(Polytope{2},model),Γface_to_bgface)

order  = 0
degree = 1

reffe_nd = ReferenceFE(nedelec,Float64,order)
V  = TestFESpace(Dc2Dp3model, reffe_nd ; conformity=:HCurl)
U = TrialFESpace(V,u)
reffe = ReferenceFE(lagrangian,Float64,order)
Q = TestFESpace(Dc2Dp3model,reffe,conformity=:L2)
P = TrialFESpace(Q)
uh = FEFunction(V,rand(num_free_dofs(V)))
vh = interpolate_everywhere(uh,V)

Ω = Triangulation(Dc2Dp3model)
dΩ = Measure(Ω,2*order)
e=sqrt(sum(∫((uh-vh)⋅(uh-vh))dΩ))
@test e < 1.0e-12

# Test a case for an oriented mesh for which the unit normals on a 
# shared facet among two cells point in opposite directions from
# the perspective of either cell. Properly handling this case requires
# an implementation of Nedelec Reference FEs which incorporates global
# mesh information when evaluating the DoFs on each cell.

function setup_model() 
    ptr  = [ 1, 9, 17 ]
    data = [ 5, 1, 7, 3, 6, 2, 8, 4, 11, 9, 12, 10, 7, 3, 8, 4 ]
    cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
    node_coordinates = Vector{Point{3,Float64}}(undef,12)
    node_coordinates[1]=Point{3,Float64}(1.0,0.0,0.0)
    node_coordinates[2]=Point{3,Float64}(1.0,0.0,1.0)
    node_coordinates[3]=Point{3,Float64}(1.0,1.0,0.0)
    node_coordinates[4]=Point{3,Float64}(1.0,1.0,1.0)
    node_coordinates[5]=Point{3,Float64}(0.0,0.0,0.0)
    node_coordinates[6]=Point{3,Float64}(0.0,0.0,1.0)
    node_coordinates[7]=Point{3,Float64}(0.0,1.0,0.0)
    node_coordinates[8]=Point{3,Float64}(0.0,1.0,1.0)
    node_coordinates[9] =Point{3,Float64}(1.0,2.0,0.0)
    node_coordinates[10]=Point{3,Float64}(1.0,2.0,1.0)
    node_coordinates[11]=Point{3,Float64}(0.0,2.0,0.0)
    node_coordinates[12]=Point{3,Float64}(0.0,2.0,1.0)
    polytope=HEX
    scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
    cell_types=collect(Fill(1,length(cell_vertex_lids)))
    cell_reffes=[scalar_reffe]
    grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                            cell_vertex_lids,
                                            cell_reffes,
                                            cell_types,
                                            Gridap.Geometry.NonOriented())
    m=Gridap.Geometry.UnstructuredDiscreteModel(grid)
    m  
end

model = setup_model()
order = 1 # Smallest nedelec FE order for which we have DoFs on the cell facets
R = TestFESpace(model, ReferenceFE(nedelec,Float64,order);conformity=:Hcurl)

f(x) = VectorValue(x[2]*x[3],x[1]*x[3],x[1]*x[2])
fh = interpolate(f,R)
Ω = Triangulation(model)
dΩ = Measure(Ω,2*order+1)
eh = fh - f
@test_broken sum(∫(eh⋅eh)dΩ) < 1.0e-12

#using Gridap.Visualization
#writevtk(trian,"trian",nsubcells=10,cellfields=["uh"=>uh])


end # module
