module RaviartThomasRefFEsTest

using Test
using Gridap.Polynomials
using Gridap.Fields
using Gridap.TensorValues
using Gridap.Fields: MockField
using Gridap.ReferenceFEs

p = QUAD
D = num_dims(QUAD)
et = Float64
order = 0

reffe = RaviartThomasRefFE(et,p,order)
test_reference_fe(reffe)
@test num_terms(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 4
@test Conformity(reffe) == DivConformity()

p = QUAD
D = num_dims(QUAD)
et = Float64
order = 1

reffe = RaviartThomasRefFE(et,p,order)
test_reference_fe(reffe)
@test num_terms(get_prebasis(reffe)) == 12
@test num_dofs(reffe) == 12
@test get_order(get_prebasis(reffe)) == 1

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

# Factory function
reffe = ReferenceFE(QUAD,:RaviartThomas,0)
@test num_terms(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 4
@test Conformity(reffe) == DivConformity()

reffe = ReferenceFE(QUAD,:RaviartThomas,Float64,0)
@test num_terms(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 4
@test Conformity(reffe) == DivConformity()

@test Conformity(reffe,:L2) == L2Conformity()
@test Conformity(reffe,:Hdiv) == DivConformity()
@test Conformity(reffe,:HDiv) == DivConformity()

end # module
