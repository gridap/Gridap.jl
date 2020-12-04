module DiscontinuousFESpacesTests

using Test
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces

domain =(0,1,0,1)
partition = (3,3)
model = CartesianDiscreteModel(domain,partition)

order = 3

reffe = ReferenceFE(:Lagrangian,Float64,order)
V = FESpace(model,reffe,conformity=:L2)
@test isa(V,UnconstrainedFESpace)
@test isa(V.cell_dofs_ids,Table{Int32})
test_single_field_fe_space(V)
@test num_free_dofs(V) == num_cells(model)*(order+1)^2

U = V
f(x) = sin(pi*x[1])*cos(2*pi*x[2])
fh = interpolate(f,U)
uh = FEFunction(V,rand(num_free_dofs(V)))

reffe = ReferenceFE(:Lagrangian,Float64,order;space=:S)
V = FESpace(model,reffe,conformity=:L2)
@test isa(V,UnconstrainedFESpace)
@test isa(V.cell_dofs_ids,Table{Int32})
test_single_field_fe_space(V)

U = V
fh = interpolate(f,U)
uh = FEFunction(V,rand(num_free_dofs(V)))

end # module
