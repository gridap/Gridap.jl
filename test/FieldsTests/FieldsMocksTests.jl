module FieldsMocksTests

using Gridap

using ..FieldsMocks

f = MockField(2,Float64)

p1 = Point(1,2)
p2 = Point(0,3)
p3 = Point(7,4)
p4 = Point(8,1)
p = [p1,p2,p3,p4]
v = Float64[1,0,7,8]
g1 = VectorValue(1.0,0.0)
g = fill(g1,4)

test_field(f,p,v,g)

b = MockBasis(2,Float64)

v = zeros(Float64,3,4)
v1 = Float64[1,0,7,8]
v[1,:] = 1*v1
v[2,:] = 2*v1
v[3,:] = 3*v1

g = zeros(VectorValue{2,Float64},3,4)
g1 = fill(g1,4)
g[1,:] = 1*g1
g[2,:] = 2*g1
g[3,:] = 3*g1

test_basis(b,p,v,g)

f = MockGeomap(2,Int)
v= 3 .* p

g1 = 3*TensorValue(1,0,0,1)
g = fill(g1,4)

test_field(f,p,v,g)

end # module
