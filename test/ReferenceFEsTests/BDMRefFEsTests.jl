module BDMRefFEsTest

using Test
using Gridap.Polynomials
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Fields: MockField
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Gridap.Arrays
using Gridap.Visualization

using Gridap
using Gridap.FESpaces

p = TRI
D = num_dims(TRI)
et = Float64
order = 1

reffe = BDMRefFE(et,p,order)

@test length(get_prebasis(reffe)) == 6
@test get_order(get_prebasis(reffe)) == 1
@test num_dofs(reffe) == 6
@test Conformity(reffe) == DivConformity()

prebasis = get_prebasis(reffe)
dof_basis = get_dof_basis(reffe)

v = VectorValue(3.0,0.0)
field = GenericField(x->v*x[1])

cache = return_cache(dof_basis,field)
r = evaluate!(cache, dof_basis, field)
test_dof_array(dof_basis,field,r)

cache = return_cache(dof_basis,prebasis)
r = evaluate!(cache, dof_basis, prebasis)
test_dof_array(dof_basis,prebasis,r)


p = TET
D = num_dims(TET)
et = Float64
order = 1

reffe = BDMRefFE(et,p,order)
test_reference_fe(reffe)
@test length(get_prebasis(reffe)) == 12
@test num_dofs(reffe) == 12
@test get_order(get_prebasis(reffe)) == 1
@test Conformity(reffe) == DivConformity()

p = TET
D = num_dims(p)
et = Float64
order = 2

reffe = BDMRefFE(et,p,order)
test_reference_fe(reffe)
@test length(get_prebasis(reffe)) == 30
@test num_dofs(reffe) == 30
@test get_order(get_prebasis(reffe)) == 2
@test Conformity(reffe) == DivConformity()

prebasis = get_prebasis(reffe)
dof_basis = get_dof_basis(reffe)

v = VectorValue(0.0,3.0,0.0)
field = GenericField(x->v)

cache = return_cache(dof_basis,field)
r = evaluate!(cache, dof_basis, field)
test_dof_array(dof_basis,field,r)

cache = return_cache(dof_basis,prebasis)
r = evaluate!(cache, dof_basis, prebasis)
test_dof_array(dof_basis,prebasis,r)

# Factory function
reffe = ReferenceFE(TET,bdm,1)
@test length(get_prebasis(reffe)) == 12
@test get_order(get_prebasis(reffe)) == 1
@test num_dofs(reffe) == 12
@test Conformity(reffe) == DivConformity()

reffe = ReferenceFE(TET,bdm,Float64,1)
@test length(get_prebasis(reffe)) == 12
@test get_order(get_prebasis(reffe)) == 1
@test num_dofs(reffe) == 12
@test Conformity(reffe) == DivConformity()

@test Conformity(reffe,:L2) == L2Conformity()
@test Conformity(reffe,:Hdiv) == DivConformity()
@test Conformity(reffe,:HDiv) == DivConformity()

@test BDM() == bdm

end # module
