module BlockFieldArrays


using Gridap.Fields
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields: MaskedField, MaskedFieldArray

using Test

p = Point(1.0,2.0)
np = 4
x = fill(p,np)
z = fill(p,0)

v = VectorValue{2,Float64}[(1,1),(4,2),(3,5)]
f = MockField.(v)

fi = first(f)

mi = MaskedField(fi,false)
test_field(mi,p,fi(p))
test_field(mi,x,fi(x))

mi = MaskedField(fi,true)
test_field(mi,p,0*fi(p))
test_field(mi,x,0*fi(x))

m = MaskedFieldArray(axes(f),f,false)
test_field_array(m,p,evaluate(f,p))
test_field_array(m,x,evaluate(f,x))

m = MaskedFieldArray(axes(f),f,true)
test_field_array(m,p,0*evaluate(f,p))
test_field_array(m,x,0*evaluate(f,x))


v1 = VectorValue{2,Float64}[(1,1),(4,2),(3,5)]
f1 = MockField.(v1)
v2 = Float64[1,2,3,4]
f2 = MockField.(v2)

axs = (append_ranges([axes(f1)[1],axes(f2)[1]]),)

f11 = MaskedFieldArray(axes(f1),f1,false)
f12 = MaskedFieldArray(axes(f2),f1,true)

b1 = BlockArrayCoo(axs,[(1,)],[f11],[1,-1],[f12])

display(b1)

display(evaluate(b1,p))


kk



display(m)
#test_field_array(m,p,evaluate(f,p))
#test_field_array(m,x,evaluate(f,x))

m = MaskedFieldArray(axes(f),f,true)
#test_field_array(m,p,0*evaluate(f,p))
#test_field_array(m,x,0*evaluate(f,x))




#v1 = VectorValue{2,Float64}[(1,1),(4,2),(3,5)]
#f1 = MockField.(v1)
#
#v2 = Float64[1,2,3]
#f2 = MockField.(v1)


end # module
