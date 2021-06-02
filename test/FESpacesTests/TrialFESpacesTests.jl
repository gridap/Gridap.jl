module TrialFESpacesTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Fields
using Gridap.FESpaces
using Gridap.CellData

domain =(0,1,0,1,0,1)
partition = (3,3,3)
model = CartesianDiscreteModel(domain,partition)

order = 2
reffe = ReferenceFE(lagrangian,Float64,order)
V = FESpace(model,reffe,dirichlet_tags=["tag_01","tag_10"])

U = TrialFESpace(V,[4,3])
@test get_cell_is_dirichlet(U) === U.space.cell_is_dirichlet
@test U.dirichlet_values == compute_dirichlet_values_for_tags(V,[4,3])
v = copy(U.dirichlet_values)
v .= 0
isa(V,SingleFieldFESpace)
ud = compute_dirichlet_values_for_tags!(v,copy(v),V,[4,3])
@test all(ud .== v)
test_single_field_fe_space(U)
U = TrialFESpace!(v,V,[4,3])

matvecdata = ([],[],[])
matdata = ([],[],[])
vecdata = ([],[])
test_single_field_fe_space(U,matvecdata,matdata,vecdata)

@test get_dirichlet_dof_values(U) == [4.0, 3.0, 3.0, 3.0, 3.0, 3.0]
TrialFESpace!(U,[1,2])
@test get_dirichlet_dof_values(U) == [1.0, 2.0, 2.0, 2.0, 2.0, 2.0]

U0 = HomogeneousTrialFESpace(V)
@test get_dirichlet_dof_values(U0) == zeros(6)

U0 = HomogeneousTrialFESpace!(v,V)
@test v === get_dirichlet_dof_values(U0)
@test v == zeros(6)
@test get_dirichlet_dof_values(U0) == zeros(6)

u(x) = x[1]
U = TrialFESpace(V,u)
uh = interpolate(u,U)
e = u - uh
trian = Triangulation(model)
quad = CellQuadrature(trian,order)

el2 = sqrt(sum(integrate(inner(e,e),quad)))
@test el2 < 1.0e-10

uh = zero(U)
cellidsL = [4,2,1,3]
cellidsR = [2,4,3,1]
cellidsS = SkeletonPair(cellidsL,cellidsR)
cell_vals = get_cell_dof_values(uh,cellidsS)
@test isa(cell_vals[1],ArrayBlock)

cell_dofs = get_cell_dof_ids(U,cellidsS)
@test isa(cell_dofs[1],ArrayBlock)

U0 = HomogeneousTrialFESpace(U)
@test get_dirichlet_dof_values(U0) == zeros(6)

#trian = get_triangulation(model)
#
#using Gridap.Visualization
#
#writevtk(trian,"trian",cellfields=["uh"=>uh])

end # module
