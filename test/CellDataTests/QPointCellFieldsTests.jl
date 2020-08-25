module QPointCellFieldsTests

using Test
using FillArrays
using Gridap.Helpers
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
ϕl = lincomb(fl,cl)

degree = 3
quad = CellQuadrature(degree,[QUAD],Fill(1,l))
q = get_coordinates(quad)

s_q = [i*ones(size(qi)) for (i,qi) in enumerate(q)]
a = CellData.ArrayOfEvaluatedFields(s_q,q)
test_array_of_fields(a,q,s_q)

v = VectorValue(3.0,4.0)
cf = QPointCellField(v,ϕl,quad)
cf_q = evaluate(cf,q)
test_cell_field(cf,q,cf_q)

cf = CellField(v,ϕl,quad)
cf_q = evaluate(cf,q)
test_cell_field(cf,q,cf_q)

cf = CellField(1.0,ϕl,quad)
@test sum(integrate(cf,ϕl,quad)) ≈ 213.33333333333323

a = CellField(1.0,ϕl,quad)
b = CellField(1.0,ϕl,quad)
c = CellField(1.0,ϕl,quad)

function updater(a,b,c)
  b_new = 2*b
  c_new = 4*c
  true, b_new, c_new
end

@test evaluate(a,q)[1] == fill(1.0,4)
@test evaluate(b,q)[1] == fill(1.0,4)
@test evaluate(c,q)[1] == fill(1.0,4)

update_state_variables!(updater,quad,a,b,c)
@test evaluate(a,q)[1] == fill(1.0,4)
@test evaluate(b,q)[1] == fill(2.0,4)
@test evaluate(c,q)[1] == fill(4.0,4)

update_state_variables!(updater,quad,a,b,c)
@test evaluate(a,q)[1] == fill(1.0,4)
@test evaluate(b,q)[1] == fill(4.0,4)
@test evaluate(c,q)[1] == fill(16.0,4)

end # module
