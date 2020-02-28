module ExtendedFESpacesTests

using Test
using Gridap.Arrays
using Gridap.Algebra
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.Integration
using Gridap.Fields

n = 10
mesh = (n,n)
domain = (0,1,0,1)
order = 1
model = CartesianDiscreteModel(domain, mesh)

trian = Triangulation(model)

const R = 0.7

function is_in(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  d = x[1]^2 + x[2]^2 - R^2
  d < 0
end

oldcell_to_coods = get_cell_coordinates(trian)

incell_to_cell = findall(collect1d(apply(is_in,oldcell_to_coods)))

trian_in = RestrictedTriangulation(trian, incell_to_cell)

reffes = [LagrangianRefFE(Float64,get_polytope(p),order) for p in get_reffes(trian_in)]

V_in = DiscontinuousFESpace(reffes,trian_in)

V = ExtendedFESpace(V_in, trian_in)
test_single_field_fe_space(V)

U = TrialFESpace(V)
test_single_field_fe_space(U)

u(x) = x[1]+x[2]

uh = interpolate(U,u)

uh_in = restrict(uh,trian_in)


using Gridap.Visualization
writevtk(trian_in,"trian_in",cellfields=["uh"=>uh_in])
writevtk(trian,"trian")#,cellfields=["uh"=>uh])




end # module
