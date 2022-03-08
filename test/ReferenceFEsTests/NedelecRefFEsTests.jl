module NedelecRefFEsTest

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

reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test num_terms(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 4

@test Conformity(reffe) == CurlConformity()

p = QUAD
D = num_dims(QUAD)
et = Float64
order = 1

reffe = NedelecRefFE(et,p,order)
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

p = TET
D = num_dims(p)
et = Float64
order = 0

reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test num_terms(get_prebasis(reffe)) == 6
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 6
@test Conformity(reffe) == CurlConformity()

p = TRI
D = num_dims(p)
et = Float64
order = 0

reffe = NedelecRefFE(et,p,order)
test_reference_fe(reffe)
@test num_terms(get_prebasis(reffe)) == 3
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 3
@test Conformity(reffe) == CurlConformity()
dof_basis = get_dof_basis(reffe)
display(dof_basis.nodes)
display(dof_basis.face_nodes)
display(dof_basis.face_moments)
display(get_face_own_dofs(reffe))
display(get_dimranges(get_polytope(reffe)))

#using Gridap.Geometry
#using Gridap.Visualization
#grid = compute_reference_grid(p,10)
#x = get_node_coordinates(grid)
#shapes = get_shapefuns(reffe)
#ux = evaluate(shapes,x)
#writevtk(grid,"nede",nodaldata=["s$i"=>ux[:,i] for i in 1:num_dofs(reffe)])


# Factory function
reffe = ReferenceFE(QUAD,nedelec,0)
@test num_terms(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 4
@test Conformity(reffe) == CurlConformity()

reffe = ReferenceFE(QUAD,nedelec,Float64,0)
@test num_terms(get_prebasis(reffe)) == 4
@test get_order(get_prebasis(reffe)) == 0
@test num_dofs(reffe) == 4
@test Conformity(reffe) == CurlConformity()

@test Conformity(reffe,:L2) == L2Conformity()
@test Conformity(reffe,:Hcurl) == CurlConformity()
@test Conformity(reffe,:HCurl) == CurlConformity()

@test Nedelec() == nedelec

end # module
