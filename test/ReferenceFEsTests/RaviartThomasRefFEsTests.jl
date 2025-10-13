module RaviartThomasRefFEsTest

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

p = QUAD
D = num_dims(QUAD)
et = Float64
order = 0

reffe = RaviartThomasRefFE(et,p,order)
test_reference_fe(reffe)
@test length(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 1
@test num_dofs(reffe) == 4
@test Conformity(reffe) == DivConformity()

@test_warn "falling back to `sh_is_pb=false`" RaviartThomasRefFE(et,p,order; sh_is_pb=true, poly_type=Monomial)

p = QUAD
D = num_dims(QUAD)
et = Float64
order = 1

reffe = RaviartThomasRefFE(et,p,order)
test_reference_fe(reffe)
@test length(get_prebasis(reffe)) == 12
@test num_dofs(reffe) == 12
@test get_order(get_prebasis(reffe)) == 2

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

###

order = 0
p = TRI
D = num_dims(TRI)
et = Float64

reffe = RaviartThomasRefFE(et,p,order)
dof_basis = get_dof_basis(reffe)
prebasis = get_prebasis(reffe)

# By default, on simplices, the element uses a raw poly basis as shapefuns, and a dof prebasis
predofs = get_dof_basis(reffe).predofs
nodes, nf_nodes, nf_moments =  get_nodes(predofs), get_face_nodes_dofs(predofs), get_face_moments(predofs)

order = 3
p = TRI
D = num_dims(TRI)
et = Float64

reffe = RaviartThomasRefFE(et,p,order)
prebasis = get_prebasis(reffe)
dof_basis = get_dof_basis(reffe)

v = VectorValue(3.0,0.0)
field = GenericField(x->v*x[1])

predofs = get_dof_basis(reffe).predofs
nodes, nf_nodes, nf_moments =  get_nodes(predofs), get_face_nodes_dofs(predofs), get_face_moments(predofs)

cache = return_cache(dof_basis,field)
r = evaluate!(cache, dof_basis, field)
test_dof_array(dof_basis,field,r)

cache = return_cache(dof_basis,prebasis)
r = evaluate!(cache, dof_basis, prebasis)
test_dof_array(dof_basis,prebasis,r)

###

p = TET
D = num_dims(TET)
et = Float64
order = 0

reffe = RaviartThomasRefFE(et,p,order)
test_reference_fe(reffe)
@test length(get_prebasis(reffe)) == 4
@test num_dofs(reffe) == 4
@test get_order(get_prebasis(reffe)) == 1
@test Conformity(reffe) == DivConformity()


p = TET
D = num_dims(p)
et = Float64
order = 2

reffe = RaviartThomasRefFE(et,p,order)
test_reference_fe(reffe)
@test length(get_prebasis(reffe)) == 36
@test num_dofs(reffe) == 36
@test get_order(get_prebasis(reffe)) == 3
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
reffe = ReferenceFE(QUAD,raviart_thomas,0)
@test reffe == ReferenceFE(QUAD,:Q⁻,1,1; rotate_90=true) # r=1, k=1
@test length(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 1
@test num_dofs(reffe) == 4
@test Conformity(reffe) == DivConformity()

@test_warn "falling back to `sh_is_pb=false`" ReferenceFE(TET,raviart_thomas,0; poly_type=Monomial)

reffe = ReferenceFE(QUAD,raviart_thomas,Float64,0)
@test reffe == ReferenceFE(QUAD,:Q⁻,1,1, Float64; rotate_90=true) # r=1, k=1
@test length(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 1
@test num_dofs(reffe) == 4
@test Conformity(reffe) == DivConformity()

reffe = ReferenceFE(HEX,raviart_thomas,0)
@test reffe == ReferenceFE(HEX,:Q⁻,1,2) # r=1, k=2

reffe = ReferenceFE(TET,raviart_thomas,0)
@test reffe == ReferenceFE(TET,:P⁻,1,2) # r=1, k=2

@test Conformity(reffe,:L2) == L2Conformity()
@test Conformity(reffe,:Hdiv) == DivConformity()
@test Conformity(reffe,:HDiv) == DivConformity()

@test RaviartThomas() == raviart_thomas

end # module
