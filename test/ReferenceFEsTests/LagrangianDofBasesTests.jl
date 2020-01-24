module LagrangianDofBasesTests

using Test

using Gridap.TensorValues
using Gridap.Fields
using Gridap.Fields: MockField, MockBasis
using Gridap.Fields: MockField
using Gridap.ReferenceFEs

x1 = Point(0,0)
x2 = Point(1,0)
x3 = Point(0,1)
x4 = Point(1,1)
x = [x1,x2,x3,x4]
np = length(x)

db = LagrangianDofBasis(Float64,x)

@test db.nodes === x
@test db.node_and_comp_to_dof == collect(1:np)
@test db.dof_to_node == collect(1:np)
@test db.dof_to_comp == fill(1,np)

v = 3.0
d = 2
f = MockField{d}(v)
fx = evaluate(f,x)
test_dof(db,f,fx)

ndof = 8
b = MockBasis{d}(v,ndof)
bx = evaluate(b,x)
test_dof(db,b,bx)

T = VectorValue{2,Float64}
db = LagrangianDofBasis(T,x)
@test db.nodes === x
@test db.node_and_comp_to_dof == VectorValue{2,Int}[(1, 5), (2, 6), (3, 7), (4, 8)]
@test db.dof_to_node == [1, 2, 3, 4, 1, 2, 3, 4]
@test db.dof_to_comp == [1, 1, 1, 1, 2, 2, 2, 2]

v = VectorValue(1,2)
d = 2
f = MockField{d}(v)
dbf = [0, 1, 0, 1, 0, 2, 0, 2]
test_dof(db,f,dbf)

ndof = 8
b = MockBasis{d}(v,ndof)
bx = evaluate(b,x)
dbb = [
  0 0 0 0 0 0 0 0; 1 1 1 1 1 1 1 1; 0 0 0 0 0 0 0 0; 1 1 1 1 1 1 1 1;
  0 0 0 0 0 0 0 0; 2 2 2 2 2 2 2 2; 0 0 0 0 0 0 0 0; 2 2 2 2 2 2 2 2]
test_dof(db,b,dbb)

end # module
