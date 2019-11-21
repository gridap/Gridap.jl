module FieldInterfaceTests

using Test
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Fields: MockField

v = 3
w = 4.0
d = 2
xi = Point(3,4)
np = 4
x = fill(xi,np)
f = MockField{d}(v)
g = MockField{d}(w)

cf = field_cache(f,x)
cg = field_cache(g,x)
c = field_caches((f,g),x)
@test (cf,cg) == c

fx = evaluate_field!(cf,f,x)
gx = evaluate_field!(cg,g,x)
fgx = evaluate_fields!(c,(f,g),x)
@test (fx,gx) == fgx

Tf = field_return_type(f,x)
Tg = field_return_type(g,x)
Tfg = field_return_types((f,g),x)
@test (Tf,Tg) == Tfg

∇f = field_gradient(f)
∇g = field_gradient(g)
∇fg = field_gradients(f,g)
@test (∇f,∇g) == ∇fg

x = Point(1,2)
@test gradient_type(Float64,x) == VectorValue{2,Float64}
@test gradient_type(VectorValue{2,Float64},x) == TensorValue{2,Float64,4}
@test gradient_type(VectorValue{3,Float64},x) == MultiValue{Tuple{2,3},Float64,2,6}

end # module
