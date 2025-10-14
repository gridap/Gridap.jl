module ModalScalarRefFEsTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.Io

nodal = false

reffe = ModalScalarRefFE(Float64,QUAD,2; F=:S)
@test reffe == ReferenceFE(QUAD,:S,2,0)
@test reffe == ReferenceFE(QUAD,:S,2,0; nodal)
test_reference_fe(reffe)

reffe = ModalScalarRefFE(Float64,QUAD,4; F=:S)
@test reffe == ReferenceFE(QUAD,:S,4,0)
@test reffe == ReferenceFE(QUAD,:S,4,0; nodal)
test_reference_fe(reffe)

reffe = ModalScalarRefFE(Float64,HEX,2; F=:S)
@test reffe == ReferenceFE(HEX,:S,2,0)
@test reffe == ReferenceFE(HEX,:S,2,0; nodal)
test_reference_fe(reffe)

p = get_polytope(reffe)
test_polytope(p)
@test LagrangianRefFE(p) == HEX8

poly_type=Polynomials.ModalC0
sh_is_pb=true

order = 4
reffe = ModalScalarRefFE(Float64,HEX,order; F=:S)
test_reference_fe(reffe)
@test reffe == ReferenceFE(HEX,:S,4,0; nodal)
@test get_order(reffe) == order

@test_throws "hierarchical" reffe = ModalScalarRefFE(Float64,HEX,order; F=:S, poly_type=Bernstein)

reffe = ModalScalarRefFE(Float64,HEX,order; F=:S, poly_type)
test_reference_fe(reffe)
@test reffe == ReferenceFE(HEX,:S,4,0; poly_type)
@test reffe == ReferenceFE(HEX,:S,4,0; poly_type, nodal)

reffe = ModalScalarRefFE(Float64,HEX,order; F=:S, poly_type, sh_is_pb)
test_reference_fe(reffe)
@test reffe == ReferenceFE(HEX,:S,4,0; poly_type, sh_is_pb)
@test reffe == ReferenceFE(HEX,:S,4,0; poly_type, sh_is_pb, nodal)

@test_warn "falling back to `sh_is_pb=false`" ModalScalarRefFE(Float64,HEX,order; F=:S, poly_type=Monomial, sh_is_pb)

end # module
