module FESpaceFactoriesTests

using Test

using Gridap.Geometry
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.CellData

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

reffe = ReferenceFE(lagrangian,Float64,order)

V = FESpace(model,reffe)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,reffe,conformity=:H1)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,reffe,conformity=:L2)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,reffe,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

# From a cell-wise vector of reffes

cell_reffe = ReferenceFE(model,lagrangian,Float64,order)

V = FESpace(model,cell_reffe)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,cell_reffe,conformity=:H1)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,cell_reffe,conformity=:L2)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,cell_reffe,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

# From a CellFE

cell_fe = FiniteElements(PhysicalDomain(),model,lagrangian,Float64,order)

V = FESpace(model,cell_fe)
@test isa(V,UnconstrainedFESpace)


cell_fe = FiniteElements(PhysicalDomain(),model,lagrangian,Float64,order,conformity=:L2)
V = FESpace(model,cell_fe)
@test isa(V,UnconstrainedFESpace)

V = FESpace(model,cell_fe,constraint=:zeromean)
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
Ω_in = Triangulation(model,cell_mask)

# From a single reffe

V = FESpace(Ω_in,QUAD4)
@test isa(V,UnconstrainedFESpace)

V = FESpace(Ω_in,QUAD4,conformity=:L2)
@test isa(V,UnconstrainedFESpace)

V = FESpace(Ω_in,QUAD4,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

# From parameter list describing the reffe

reffe = ReferenceFE(lagrangian,Float64,order)

V = FESpace(Ω_in,reffe)
@test isa(V,UnconstrainedFESpace)

V = FESpace(Ω_in,reffe,conformity=:L2)
@test isa(V,UnconstrainedFESpace)

V = FESpace(Ω_in,reffe,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

# From a cell-wise vector of reffes

model_in = get_active_model(Ω_in)
cell_reffe = ReferenceFE(model_in,lagrangian,Float64,order)

V = FESpace(model_in,cell_reffe,trian=Ω_in)
@test isa(V,UnconstrainedFESpace)
@test Ω_in === get_triangulation(V)

V = FESpace(model_in,cell_reffe,trian=Ω_in,conformity=:L2)
@test isa(V,UnconstrainedFESpace)
@test Ω_in === get_triangulation(V)

V = FESpace(model_in,cell_reffe,trian=Ω_in,conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)
@test Ω_in === get_triangulation(V)

# From a CellFE

cell_fe = FiniteElements(PhysicalDomain(),model_in,lagrangian,Float64,order)

V = FESpace(model_in,cell_fe,trian=Ω_in)
@test isa(V,UnconstrainedFESpace)
@test Ω_in === get_triangulation(V)

cell_fe = FiniteElements(PhysicalDomain(),model_in,lagrangian,Float64,order,conformity=:L2)
V = FESpace(model_in,cell_fe,trian=Ω_in)
@test isa(V,UnconstrainedFESpace)
@test Ω_in === get_triangulation(V)

V = FESpace(model_in,cell_fe,trian=Ω_in,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)
@test Ω_in === get_triangulation(V)

Γ = Boundary(model)
reffe = ReferenceFE(lagrangian,Float64,1)
V = FESpace(Γ,reffe)
@test get_triangulation(V) === Γ
V = FESpace(Γ,reffe,conformity=:H1)
@test get_triangulation(V) === Γ
V = FESpace(Γ,reffe,conformity=:L2)
@test get_triangulation(V) === Γ

end # module
