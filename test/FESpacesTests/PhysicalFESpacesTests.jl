module PhysicalFESpacesTests

using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.Fields
using Gridap.FESpaces
using Gridap.Polynomials
using Gridap.CellData
using Gridap.TensorValues
using Test

u(x) = x[1]+2*x[2] + 3

domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
order = 1

Ω = Triangulation(model)
dΩ = Measure(Ω,order)

cell_fe = FiniteElements(ReferenceDomain(),model,lagrangian,Float64,order)
V = FESpace(model,cell_fe)
@test num_free_dofs(V) == 9
cell_fe = FiniteElements(ReferenceDomain(),model,lagrangian,Float64,order,conformity=:L2)
V = FESpace(model,cell_fe)
@test num_free_dofs(V) == 16

@test DomainStyle(get_fe_basis(V)) == ReferenceDomain()
@test DomainStyle(get_fe_dof_basis(V)) == ReferenceDomain()

uh = interpolate(u,V)
e = u - uh
@test sqrt(sum(∫(e*e)*dΩ)) < 10e-8

cell_fe = FiniteElements(PhysicalDomain(),model,lagrangian,Float64,order)
V = FESpace(model,cell_fe)
@test num_free_dofs(V) == 9

cell_fe = FiniteElements(PhysicalDomain(),model,lagrangian,Float64,order,conformity=:L2)
V = FESpace(model,cell_fe)
@test num_free_dofs(V) == 16

@test DomainStyle(get_fe_basis(V)) == PhysicalDomain()
@test DomainStyle(get_fe_dof_basis(V)) == PhysicalDomain()

uh = interpolate(u,V)
e = u - uh
@test sqrt(sum(∫(e*e)*dΩ)) < 10e-8

end #module
