module FESpaceFactoriesTests

using Test

using Gridap.Geometry
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.CellData
using Gridap.FESpaces: ExtendedFESpace

domain = (0,1,0,1)
partition = (10,10)
model = CartesianDiscreteModel(domain,partition)
order = 3

# From a single reffe

V = FESpace(model,QUAD4)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,QUAD4,conformity=:H1)
@test isa(V,UnconstrainedFESpace)
@test num_free_dofs(V) == 121

V = FESpace(model,QUAD4,conformity=:L2)
@test isa(V,UnconstrainedFESpace)
@test num_free_dofs(V) == 400

V = FESpace(model,QUAD4,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

# From parameter list describing the reffe

reffe = ReferenceFE(:Lagrangian,Float64,order)

V = FESpace(model,reffe)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,reffe,conformity=:H1)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,reffe,conformity=:L2)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,reffe,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

# From a cell-wise vector of reffes

cell_reffe = ReferenceFE(model,:Lagrangian,Float64,order)

V = FESpace(model,cell_reffe)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,cell_reffe,conformity=:H1)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,cell_reffe,conformity=:L2)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,cell_reffe,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

# From a CellFE

cell_fe = FiniteElements(PhysicalDomain(),model,:Lagrangian,Float64,order)

V = FESpace(model,cell_fe)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,cell_fe,conformity=:L2)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,cell_fe,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

# Now from a restricted model

const R = 0.7

function is_in(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  d = x[1]^2 + x[2]^2 - R^2
  d < 0
end

Ω = Triangulation(model)
cell_coords = get_cell_coordinates(Ω)
cell_mask = lazy_map(is_in,cell_coords)
model_in = DiscreteModel(model,cell_mask)

# From a single reffe

V = FESpace(model_in,QUAD4)
@test isa(V,ExtendedFESpace)

V = FESpace(model_in,QUAD4,conformity=:L2)
@test isa(V,ExtendedFESpace)

V = FESpace(model_in,QUAD4,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

# From parameter list describing the reffe

reffe = ReferenceFE(:Lagrangian,Float64,order)

V = FESpace(model_in,reffe)
@test isa(V,ExtendedFESpace)

V = FESpace(model_in,reffe,conformity=:L2)
@test isa(V,ExtendedFESpace)

V = FESpace(model_in,reffe,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

# From a cell-wise vector of reffes

cell_reffe = ReferenceFE(model_in,:Lagrangian,Float64,order)

V = FESpace(model_in,cell_reffe)
@test isa(V,ExtendedFESpace)

V = FESpace(model_in,cell_reffe,conformity=:L2)
@test isa(V,ExtendedFESpace)

V = FESpace(model_in,cell_reffe,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

# From a CellFE

cell_fe = FiniteElements(PhysicalDomain(),model_in,:Lagrangian,Float64,order)

V = FESpace(model_in,cell_fe)
@test isa(V,ExtendedFESpace)

V = FESpace(model_in,cell_fe,conformity=:L2)
@test isa(V,ExtendedFESpace)

V = FESpace(model_in,cell_fe,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

end # module
