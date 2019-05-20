
include("FESpaces.jl")

module FESpacesTests

using Test

using Gridap
using Gridap.RefFEs
using Gridap.Polytopes
using Gridap.Geometry
using Gridap.Geometry.Cartesian

using ..FESpaces

model = CartesianDiscreteModel(domain=(0.0,1.0,-1.0,2.0), partition=(2,2))

D = pointdim(model)

grid = Grid(model,D)
trian = triangulation(grid)
graph = FullGridGraph(model)
labels = FaceLabels(model)
tags = [1,2,3,4]

order = 1
orders = fill(order,D)
polytope = Polytope(fill(HEX_AXIS,D)...)
fe = LagrangianRefFE{D,Float64}(polytope, orders)

fespace = ConformingFESpace(fe,trian,graph,labels,tags)

@test num_free_dofs(fespace) == 5
@test num_fixed_dofs(fespace) == 4

r = [[-1, 1, 2, 3], [1, -2, 3, 4], [2, 3, -3, 5], [3, 4, 5, -4]]

@test r == collect(fespace.cell_eqclass)

order = 2
orders = fill(order,D)
polytope = Polytope(fill(HEX_AXIS,D)...)
fe = LagrangianRefFE{D,Float64}(polytope, orders)

tags = [1,2,3,4,6,5]
fespace = ConformingFESpace(fe,trian,graph,labels,tags)

@test num_free_dofs(fespace) == 15
@test num_fixed_dofs(fespace) == 10

r = [[-1, -2, 1, 2, -7, 4, 5, 6, 12], [-2, -3, 2, 3, -8, 7, 6, 8, 13],
     [1, 2, -4, -5, 4, -9, 9, 10, 14], [2, 3, -5, -6, 7, -10, 10, 11, 15]]

@test r == collect(fespace.cell_eqclass)

end # module FESpacesTests
