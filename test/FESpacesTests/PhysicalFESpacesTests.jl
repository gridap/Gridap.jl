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

# domain = (0,1)
# partition = (3,)
domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
order = 1

trian = get_triangulation(model)
quad = CellQuadrature(trian,order)

u(x) = x[1]
T = Float64

u(x) = VectorValue(x[1],x[2])
T = VectorValue{2,Float64}

# Juno.@run TestFESpace( reffe = :Lagrangian, conformity = :H1, valuetype=T, order = order, model = model, dirichlet_tags = [1,6], dof_space = :physical)
# Juno.@enter TestFESpace( reffe = :Lagrangian, conformity = :H1, valuetype=T, order = order, model = model, dirichlet_tags = [1,6], dof_space = :physical)
# Juno.@enter TestFESpace( reffe = :Lagrangian, conformity = :H1, valuetype=T, order = order, model = model, dirichlet_tags = [1,6])

Vp = TestFESpace(
  reffe = :Lagrangian,
  conformity = :H1,
  valuetype=T,
  order = order,
  model = model,
  dirichlet_tags = [1,6],
  dof_space = :physical)

V = TestFESpace(
  reffe = :Lagrangian,
  conformity = :H1,
  valuetype=VectorValue{2,Float64},
  order = order,
  model = model,
  dirichlet_tags = [1,6])

# Vp = TestFESpace(
#   reffe = :RaviartThomas,
#   conformity = :Hdiv,
#   order = order,
#   model = model,
#   dirichlet_tags = [1,6],
#   dof_space = :physical)
#
# V = TestFESpace(
#   reffe = :RaviartThomas,
#   conformity = :Hdiv,
#   order = order,
#   model = model,
#   dirichlet_tags = [1,6])
# test_single_field_fe_space(V)

Up = TrialFESpace(Vp,u)
U = TrialFESpace(V,u)

# FEM = Gridap.FESpaces
# Juno.@enter TrialFESpace(Vp,u)
# dirichlet_values = FEM.compute_dirichlet_values_for_tags(V,u)

q = get_coordinates(quad)
r = evaluate(U.cell_basis,q)
rp = evaluate(Up.cell_basis,q)

@test all(collect(r) .≈ collect(rp))

#@test RefStyle(Up.space.cell_dof_basis) == Val{false}()
#@test RefStyle(U.space.cell_dof_basis) == Val{true}()
#
#@test RefStyle(Up.space.cell_basis) == Val{false}()
#@test RefStyle(U.space.cell_basis) == Val{true}()
#
#@test RefStyle(Vp.cell_basis) == Val{false}()
#@test RefStyle(V.cell_basis) == Val{true}()

uhf = interpolate(u,Up)
uhd = interpolate_dirichlet(u,Up)
uhdf = interpolate_everywhere(u,Up)
uhf.dirichlet_values
uhdf.dirichlet_values
uhd.dirichlet_values

uhf_r = interpolate(u,U)
uhd_r= interpolate_dirichlet(u,U)
uhdf_r = interpolate_everywhere(u,U)
uhd_r.dirichlet_values

@test uhf.free_values == uhf_r.free_values
@test uhf.dirichlet_values == uhf_r.dirichlet_values
@test uhd.dirichlet_values == uhd_r.dirichlet_values
@test uhdf.free_values == uhdf_r.free_values
@test uhdf.dirichlet_values == uhdf_r.dirichlet_values

@test uhf.free_values == uhdf.free_values
@test uhf.dirichlet_values == uhd.dirichlet_values
@test uhdf.free_values == uhf.free_values
@test uhdf.dirichlet_values == uhd.dirichlet_values
@test uhdf.dirichlet_values == uhdf_r.dirichlet_values

e = u - uhf_r

el2 = sqrt(sum(integrate(inner(e,e),quad)))

@test el2 < 1.0e-10

e = u - uhf

el2 = sqrt(sum(integrate(inner(e,e),quad)))

@test el2 < 1.0e-10

end #module
