module FESpaceWithConstantFixedTests
using Gridap
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.ReferenceFEs
using Test

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

trian = get_triangulation(model)

order = 2

V = FESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:L2)

V0 = FESpaceWithConstantFixed(V,true,rand(1:num_free_dofs(V)))
test_single_field_fe_space(V0)

@test Gridap.FESpaces.ConstantApproach(V0) == Gridap.FESpaces.FixConstant()

uh0 = interpolate(V0) do x
    sin(4*pi*(x[1]+x[2]^2)) + 3
end
using Gridap.Visualization

#writevtk(trian,"trian",nsubcells=20,cellfields=["uh0"=>uh0])

end # module
