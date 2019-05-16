module CellValuesBench

using Gridap.FieldValues
using Gridap.CellValues

include("../test/CellValuesTests/Mocks.jl")

function doloop(x)
  for xi in x
  end
end

l = 1000000

println("+++ CellValuesBench ( length = $l ) +++")

sv = 1.0
scv = ConstantCellValue(sv,l)
sv2 = 1.1
scv2 = ConstantCellValue(sv2,l)
sa = [sv, sv, sv]
sca = ConstantCellValue(sa,l)
sa2 = [sv sv; sv sv; sv sv]
sca2 = ConstantCellValue(sa2,l)

print("ConstantCellValue ->"); @time doloop(scv)
print("ConstantCellValue ->"); @time doloop(scv)

print("ConstantCellValue ->"); @time doloop(sca2)
print("ConstantCellValue ->"); @time doloop(sca2)

sv = 1.0
scv = TestCellValue(sv,l)
sv2 = 1.1
scv2 = TestCellValue(sv2,l)
sa = [sv, sv, sv]
sca = TestCellArray(sa,l)
sa2 = [sv sv; sv sv; sv sv]
sca2 = TestCellArray(sa2,l)

print("TestCellValue ->"); @time doloop(scv)
print("TestCellValue ->"); @time doloop(scv)

print("TestCellArray ->"); @time doloop(sca2)
print("TestCellArray ->"); @time doloop(sca2)

for op in (:+,:-)
  @eval begin
    scv3 = $op(scv)
    print("CellValueUnary($(string($op))) ->"); @time doloop(scv3)
    print("CellValueUnary($(string($op))) ->"); @time doloop(scv3)
  end
end

for op in (:+,:-)
  @eval begin
    scv3 = $op(sca)
    print("CellArrayUnary($(string($op))) ->"); @time doloop(scv3)
    print("CellArrayUnary($(string($op))) ->"); @time doloop(scv3)
  end
end

sca3 = cellsum(sca2,dim=2)
print("CellArrayFromCellSum ->"); @time doloop(sca3)
print("CellArrayFromCellSum ->"); @time doloop(sca3)

sca3 = cellsum(sca,dim=1)
print("CellValueFromCellSum ->"); @time doloop(sca3)
print("CellValueFromCellSum ->"); @time doloop(sca3)

sca3 = cellnewaxis(sca2,dim=2)
print("CellArrayFromNewAxis ->"); @time doloop(sca3)
print("CellArrayFromNewAxis ->"); @time doloop(sca3)

for op in (:+,:-,:*,:/,:(outer),:(inner))
  @eval begin
    scv3 = $op(scv,scv2)
    print("CellValueBinary($(string($op))) ->"); @time doloop(scv3)
    print("CellValueBinary($(string($op))) ->"); @time doloop(scv3)
  end
end

for op in (:+,:-,:*,:/,:(outer),:(inner))
  @eval begin
    sca3 = $op(sca,sca2)
    print("CellArrayBinary($(string($op))) ->"); @time doloop(sca3)
    print("CellArrayBinary($(string($op))) ->"); @time doloop(sca3)
  end
end

for op in (:+,:-,:*,:/,:(outer),:(inner))
  @eval begin
    sca3 = $op(scv,sca2)
    print("CellValueArrayBinary($(string($op))) ->"); @time doloop(sca3)
    print("CellValueArrayBinary($(string($op))) ->"); @time doloop(sca3)
  end
end

end # module CellValuesBench
