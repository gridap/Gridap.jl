module CellValuesBench

using Numa.CellValues

include("../test/CellValuesTests/Mocks.jl")

function doloop(x)
  i = 0
  for xi in x
    i += 1
  end
end

l = 1000000

println("+++ CellValuesBench ( length = $l ) +++")

sv = 1.0
scv = ConstantCellValue(sv,l)
sv2 = 1.1
scv2 = ConstantCellValue(sv2,l)
sa = [sv, sv, sv]
sca = ConstantCellArray(sa,l)
sa2 = [sv sv; sv sv; sv sv]
sca2 = ConstantCellArray(sa2,l)

print("ConstantCellValue ->"); @time doloop(scv)
print("ConstantCellValue ->"); @time doloop(scv)

print("ConstantCellArray ->"); @time doloop(sca2)
print("ConstantCellArray ->"); @time doloop(sca2)

sv = 1.0
scv = TestCellValue(sv,l)
sv2 = 1.1
scv2 = TestCellValue(sv2,l)
sa = [sv, sv, sv]
sca = TestCellValue(sa,l)
sa2 = [sv sv; sv sv; sv sv]
sca2 = TestCellValue(sa2,l)

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

for op in (:+,:-,:*,:/)
  @eval begin
    scv3 = $op(scv,scv2)
    print("CellValueBinary($(string($op))) ->"); @time doloop(scv3)
    print("CellValueBinary($(string($op))) ->"); @time doloop(scv3)
  end
end

end # module CellValuesBench
