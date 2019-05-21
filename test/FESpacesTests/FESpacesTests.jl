
include("FESpaces.jl")

module FESpacesTests

using Test

using Gridap
using Gridap.RefFEs
using Gridap.CellMaps
using Gridap.Polytopes
using Gridap.Geometry
using Gridap.Geometry.Cartesian
using Gridap.CellMaps.Testers

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
@test num_diri_dofs(fespace) == 4

r = [[-1, 1, 2, 3], [1, -2, 3, 4], [2, 3, -3, 5], [3, 4, 5, -4]]

@test r == collect(fespace.cell_eqclass)

order = 2
orders = fill(order,D)
polytope = Polytope(fill(HEX_AXIS,D)...)
fe = LagrangianRefFE{D,Float64}(polytope, orders)

tags = [1,2,3,4,6,5]
fespace = ConformingFESpace(fe,trian,graph,labels,tags)

@test num_free_dofs(fespace) == 15
@test num_diri_dofs(fespace) == 10

r = [[-1, -2, 1, 2, -7, 4, 5, 6, 12], [-2, -3, 2, 3, -8, 7, 6, 8, 13],
     [1, 2, -4, -5, 4, -9, 9, 10, 14], [2, 3, -5, -6, 7, -10, 10, 11, 15]]

@test r == collect(fespace.cell_eqclass)

fun(x) = sin(x[1])*cos(x[2])

free_vals, diri_vals = interpolated_values(fespace,fun)

rf = [0.0, 0.420735, 0.598194, 0.078012, 0.0, 0.420735,
      0.214936, 0.598194, -0.0, -0.199511, -0.283662,
      0.420735, 0.73846, -0.199511, -0.350175]

rd = [0.0, 0.259035, 0.368291, 0.420735, 0.73846, 0.151174,
      0.239713, 0.660448, 0.151174, 0.265335]

@test isapprox(free_vals,rf,rtol=1.0e-5)
@test isapprox(diri_vals,rd,rtol=1.0e-5)

uh = FEFunction(fespace,free_vals,diri_vals)

@test free_dofs(uh) === free_vals
@test diri_dofs(uh) === diri_vals
@test FESpace(uh) === fespace

cellbasis = CellBasis(fespace)

uh = interpolate(fespace,fun)
@test isa(uh,FEFunction)

quad = quadrature(trian,order=2)

q = coordinates(quad)
uhq = evaluate(uh,q)

grad_uh = gradient(uh)
grad_uhq = evaluate(grad_uh,q)

v = collect(uhq)
g = collect(grad_uhq)

test_cell_map_with_gradient(uh,q,v,g)

end # module FESpacesTests
