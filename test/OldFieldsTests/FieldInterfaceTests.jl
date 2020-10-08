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

cf = return_cache(f,x)
cg = return_cache(g,x)
c = return_caches((f,g),x)
@test (cf,cg) == c

fx = evaluate!(cf,f,x)
gx = evaluate!(cg,g,x)
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
@test return_gradient_type(Float64,x) == VectorValue{2,Float64}
@test return_gradient_type(VectorValue{2,Float64},x) == TensorValue{2,2,Float64,4}
@test return_gradient_type(VectorValue{3,Float64},x) == TensorValue{2,3,Float64,6}

end # module
