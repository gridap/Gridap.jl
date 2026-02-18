module DiscontinuousFESpacesTests

using Test
using Gridap
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

order = 3

reffe = ReferenceFE(lagrangian,Float64,order)
V = FESpace(model,reffe,conformity=:L2)
@test isa(V,UnconstrainedFESpace)
@test isa(V.cell_dofs_ids,Table{Int32})
test_single_field_fe_space(V)
@test num_free_dofs(V) == num_cells(model)*(order+1)^2

U = V
f(x) = sin(pi*x[1])*cos(2*pi*x[2])
fh = interpolate(f,U)
uh = FEFunction(V,rand(num_free_dofs(V)))

reffe = ReferenceFE(lagrangian,Float64,order;space=:S)
V = FESpace(model,reffe,conformity=:L2)
@test isa(V,UnconstrainedFESpace)
@test isa(V.cell_dofs_ids,Table{Int32})
test_single_field_fe_space(V)

U = V
fh = interpolate(f,U)
uh = FEFunction(V,rand(num_free_dofs(V)))

###################

Vp = FESpaces.PolytopalFESpace(model,Float64, 1; space=:Q)
cell_conformity = FESpaces.get_cell_conformity(Vp)
@test isa(cell_conformity, FESpaces.DiscontinuousCellConformity)
bmask = FESpaces.generate_dof_mask(Vp,get_face_labeling(model),"boundary")
@test sum(bmask) == 0
bmask_rev = FESpaces.generate_dof_mask(Vp,get_face_labeling(model),"boundary",reverse=true)
@test sum(bmask_rev) == num_free_dofs(Vp)
@test all(bmask .== .!bmask_rev)

reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},2)
V = FESpace(Triangulation(ReferenceFE{1},model),reffe;conformity=:L2,dirichlet_tags=["boundary"])
@test isa(FESpaces.get_cell_conformity(V), FESpaces.CompressedCellConformity) # I think we should change this
cell_to_tag = collect(Int8,Geometry.get_face_mask(get_face_labeling(model),"boundary",1))
cell_dofs, nfree, ndir, dirichlet_dof_tag, dirichlet_cells = FESpaces.compute_discontinuous_cell_dofs(
  FESpaces.get_cell_conformity(V), cell_to_tag, nothing
)
@test get_cell_dof_ids(V) == cell_dofs

end # module
