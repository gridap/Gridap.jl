module QPointCellFieldsTests

using Test
using FillArrays
using Gridap.Helpers
using Gridap.Arrays
using Gridap.Fields
using Gridap.CellData
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Integration

using Gridap.Fields: MockField, MockBasis, OtherMockBasis

ndofs = 4
c = fill(1.0,ndofs)
f = OtherMockBasis{2}(ndofs)

l = 10
cl = fill(c,l)
fl = Fill(f,l)
ϕ = GenericCellMap(lincomb(fl,cl))

degree = 3
ref_quad = CellQuadrature(degree,[QUAD],Fill(1,l))
q = get_coordinates(ref_quad)

s_q = [i*ones(size(qi)) for (i,qi) in enumerate(get_array(q))]
a = CellData.ArrayOfEvaluatedFields(s_q,get_array(q))
test_array_of_fields(a,get_array(q),s_q)

quad = ϕ(ref_quad)
x = get_coordinates(quad)

v = VectorValue(3.0,4.0)
cf = QPointCellField(v,quad)
cf_x = evaluate(cf,x)
test_cell_field(cf,x,cf_x)

cf = CellField(v,quad)
cf_q = evaluate(cf,x)
test_cell_field(cf,x,cf_x)

cf = CellField(1.0,quad)
@test sum(integrate(cf,quad)) ≈ 213.33333333333323
@test sum(∫(cf)*quad) ≈ 213.33333333333323

a = CellField(1.0,quad)
b = CellField(1.0,quad)
c = CellField(1.0,quad)

function updater(a,b,c)
  b_new = 2*b
  c_new = 4*c
  true, b_new, c_new
end

@test evaluate(a,x)[1] == fill(1.0,4)
@test evaluate(b,x)[1] == fill(1.0,4)
@test evaluate(c,x)[1] == fill(1.0,4)

update_state_variables!(updater,quad,a,b,c)
@test evaluate(a,x)[1] == fill(1.0,4)
@test evaluate(b,x)[1] == fill(2.0,4)
@test evaluate(c,x)[1] == fill(4.0,4)

update_state_variables!(updater,quad,a,b,c)
@test evaluate(a,x)[1] == fill(1.0,4)
@test evaluate(b,x)[1] == fill(4.0,4)
@test evaluate(c,x)[1] == fill(16.0,4)

end # module
