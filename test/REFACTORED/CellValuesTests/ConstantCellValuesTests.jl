module ConstantCellValuesTests

using Test
using Gridap
using TensorValues
using ..MapsMocks

# Constructors

l = 10
sv = 1.0
scv = ConstantCellNumber(sv,l)

sv2 = 1.1
scv2 = ConstantCellNumber(sv2,l)

test_index_cell_value( scv, fill(sv,l) )

sa = [sv, sv, sv]
sca = ConstantCellArray(sa,l)

test_index_cell_array( sca, fill(sa,l) )

sca = ConstantCellVector(sa,l)

test_index_cell_array( sca, fill(sa,l) )

sa2 = [sv sv; sv sv; sv sv]
sca2 = ConstantCellArray(sa2,l)

test_index_cell_array( sca2, fill(sa2,l) )

sca2 = ConstantCellMatrix(sa2,l)

test_index_cell_array( sca2, fill(sa2,l) )

a = VectorValue(10,10)
b = VectorValue(15,20)
p1 = VectorValue(1,1)
p2 = VectorValue(2,2)
p3 = VectorValue(3,3)
p = [p1,p2,p3]
m = MockMap(a)
r = evaluate(m,p)

cm = ConstantCellMap(m,l)
cp = ConstantCellVector(p,l)
rm = [ r for i in 1:l]
test_index_cell_map_with_index_arg(cm,cp,rm)
cr = evaluate(cm,cp)
@test isa(cr,ConstantCellArray)

scv3 = apply(-,scv,scv2)
@test isa(scv3,ConstantCellValue{Float64})
test_index_cell_value( scv3, fill(-(sv,sv2),l) )

for op in (:+,:-)
  @eval begin
    scv3 = $op(scv,scv2)
    @test isa(scv3,ConstantCellValue{Float64})
    test_index_cell_value( scv3, fill($op(sv,sv2),l) )
  end
end

sca3 = apply(-,sca2,broadcast=true)
@test isa(sca3,ConstantCellArray{Float64,2})
test_index_cell_array( sca3, fill(broadcast(-,sa2),l) )

sca3 = apply(-,sca,sca2,broadcast=true)
@test isa(sca3,ConstantCellArray{Float64,2})
test_index_cell_array( sca3, fill(broadcast(-,sa,sa2),l) )

cm = ConstantCellMap(m,l)
cp = ConstantCellVector(p,l)
rm = [ -r for i in 1:l]
cm2 = -cm
test_index_cell_map_with_index_arg(cm2,cp,rm)
cr = evaluate(cm2,cp)
@test isa(cr,ConstantCellArray)

cm = ConstantCellMap(m,l)
cp = ConstantCellVector(p,l)
rm = [ r.-p for i in 1:l]
cm2 = apply(-,cm,cp,broadcast=true)
test_index_cell_map_with_index_arg(cm2,cp,rm)
cr = evaluate(cm2,cp)
@test isa(cr,ConstantCellArray)

@test scv == scv
@test scv ≈ scv

@test cp == cp
@test cp ≈ cp

ca = ConstantCellNumber(3,5)
sv = 1.0
scv = ConstantCellNumber(sv,l)
cc = reindex(scv,ca)
r = fill(sv,5)
@assert isa(cc,ConstantCellValue)
test_iter_cell_value(cc,r)


end # module ConstantCellValuesTests
