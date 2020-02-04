module SubVectorsTests

using Gridap.Arrays

a = collect(1:10)
pini = 1
pend = 5
b = SubVector(a,pini,pend)
test_array(b,collect(pini:pend))

pini = 7
pend = 9
b = SubVector(a,pini,pend)
test_array(b,collect(pini:pend))

end # module
