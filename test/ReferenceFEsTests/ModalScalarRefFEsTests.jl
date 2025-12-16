module ModalScalarRefFEsTests

using Test
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Polynomials
using Gridap.TensorValues: SymTensorValue
using Gridap.Io

nodal = false

p, F = (TRI, :P⁻)
reffe = ModalScalarRefFE(Float64,p,2; F)
@test reffe == ReferenceFE(p,F,2,0)
@test reffe == ReferenceFE(p,F,2,0; nodal)

face_own_dofs = Vector{Int}[[1],[2],[3],[4],[5],[6],[]]
face_dofs = Vector{Int}[[1],[2],[3],[1,2,4],[1,3,5],[2,3,6],[1,2,3,4,5,6]]
@test get_face_own_dofs(reffe) == face_own_dofs
@test get_face_dofs(reffe) == face_dofs

test_reference_fe(reffe)

for (p,F) in [
  (TRI, :P⁻), (TRI, :P), (QUAD,:Q⁻), (QUAD,:S),
  (TET, :P⁻), (TET, :P), (HEX, :Q⁻), (HEX, :S),
  ]

  reffe = ModalScalarRefFE(Float64,p,1; F)
  @test reffe == ReferenceFE(p,F,1,0)
  @test reffe == ReferenceFE(p,F,1,0; nodal)
  test_reference_fe(reffe)

  reffe = ModalScalarRefFE(Float64,p,4; F)
  @test reffe == ReferenceFE(p,F,4,0)
  @test reffe == ReferenceFE(p,F,4,0; nodal)
  test_reference_fe(reffe)

  V = SymTensorValue{2,Float64}
  reffe = ModalScalarRefFE(V,p,2; F)
  @test reffe == ReferenceFE(p,F,2,0,V)
  @test reffe == ReferenceFE(p,F,2,0,V; nodal)
  test_reference_fe(reffe)
end

poly_type=Polynomials.ModalC0
change_dof=true

order = 4
reffe = ModalScalarRefFE(Float64,HEX,order; F=:S)
test_reference_fe(reffe)
@test reffe == ReferenceFE(HEX,:S,4,0; nodal)
@test get_order(reffe) == order
@test Conformity(reffe,:H1) == GradConformity()
@test Conformity(reffe,:L2) == L2Conformity()

@test_throws "hierarchical" reffe = ModalScalarRefFE(Float64,HEX,order; F=:S, poly_type=Bernstein)

reffe = ModalScalarRefFE(Float64,HEX,order; F=:S, poly_type)
test_reference_fe(reffe)
@test reffe == ReferenceFE(HEX,:S,4,0; poly_type)
@test reffe == ReferenceFE(HEX,:S,4,0; poly_type, nodal)

reffe = ModalScalarRefFE(Float64,HEX,order; F=:S, poly_type, change_dof)
test_reference_fe(reffe)
@test reffe == ReferenceFE(HEX,:S,4,0; poly_type, change_dof)
@test reffe == ReferenceFE(HEX,:S,4,0; poly_type, change_dof, nodal)

@test_warn "falling back to `change_dof=false`" ModalScalarRefFE(Float64,HEX,order; F=:S, poly_type=Monomial, change_dof)

end # module
