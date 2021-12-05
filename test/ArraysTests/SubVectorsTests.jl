module SubVectorsTests
using Gridap
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

inds = Int32[1,4,3]
a = rand(10)
sa = Gridap.Arrays.SubVector(a,1,5)
vsa = view(sa,inds)

b = rand(10)
vb = view(b,inds)

x = Base.broadcasted(+,vb,vsa)
Base.materialize!(vb,x)

end # module
