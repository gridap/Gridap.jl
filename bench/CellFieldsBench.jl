module CellFieldsBench

using Numa.FieldValues
using Numa.CellArrays
using Numa.CellFields

l = 1000000

include("../test/CellFieldsTestsMocks.jl")

function doloop(x)
  for xi in x
  end
end

println("+++ CellFieldsBench ( length = $l ) +++")

tiff = inner(tbv,tbv)

print("TensorBasisBasisInner ->"); @time doloop(tiff)
print("TensorBasisBasisInner ->"); @time doloop(tiff)

vexpand = expand(vbv,sfv2)

print("VectorScalarExpand ->"); @time doloop(vexpand)
print("VectorScalarExpand ->"); @time doloop(vexpand)

end
