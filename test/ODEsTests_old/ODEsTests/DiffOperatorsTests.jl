module DiffOperatorsTests

using Gridap
using Test
using Gridap.ODEs
using Gridap.ODEs.ODETools: ∂t, ∂tt

using ForwardDiff

f(x,t) = 5*x[1]*x[2]+x[2]^2*t^3
∂tf(x,t) = x[2]^2*3*t^2

tv = rand(Float64)
xv = Point(rand(Float64,2)...)
@test ∂tf(xv,tv) ≈ ∂t(f)(xv,tv)

F(x,t) = VectorValue([5*x[1]*x[2],x[2]^2*t^3])
∂tF(x,t) = VectorValue([0.0,x[2]^2*3*t^2])

tv = rand(Float64)
xv = Point(rand(Float64,2)...)
@test ∂tF(xv,tv) ≈ ∂t(F)(xv,tv)

# Time derivatives

f(x,t) = t^2
dtf = (x,t) -> ForwardDiff.derivative(t->f(x,t),t)
@test dtf(xv,tv) ≈ ∂t(f)(xv,tv) ≈ ∂t(f)(xv)(tv) ≈ ∂t(f)(tv)(xv)

f2(x,t) = x[1]^2
dtf2 = (x,t) -> ForwardDiff.derivative(t->f2(x,t),t)
@test dtf2(xv,tv) ≈ ∂t(f2)(xv,tv) ≈ ∂t(f2)(xv)(tv) ≈ ∂t(f2)(tv)(xv)

f2(x,t) = x[1]^t^2
dtf2 = (x,t) -> ForwardDiff.derivative(t->f2(x,t),t)
@test dtf2(xv,tv) ≈ ∂t(f2)(xv,tv) ≈ ∂t(f2)(xv)(tv) ≈ ∂t(f2)(tv)(xv)

f2(x,t) = VectorValue(x[1]^2,0.0)
dtf2 = (x,t) -> VectorValue(ForwardDiff.derivative(t -> get_array(f2(x,t)),t))
@test dtf2(xv,tv) ≈ ∂t(f2)(xv,tv) ≈ ∂t(f2)(xv)(tv) ≈ ∂t(f2)(tv)(xv)

f2(x,t) = VectorValue(x[1]^2,t)
dtf2 = (x,t) -> VectorValue(ForwardDiff.derivative(t -> get_array(f2(x,t)),t))
@test dtf2(xv,tv) ≈ ∂t(f2)(xv,tv) ≈ ∂t(f2)(xv)(tv) ≈ ∂t(f2)(tv)(xv)

f6(x,t) = TensorValue(x[1]*t,x[2]*t,x[1]*x[2],x[1]*t^2)
dtf6 = (x,t) -> TensorValue(ForwardDiff.derivative(t->get_array(f6(x,t)),t))
@test dtf6(xv,tv) ≈ ∂t(f6)(xv,tv) ≈ ∂t(f6)(xv)(tv) ≈ ∂t(f6)(tv)(xv)

# Spatial derivatives
f2(x,t) = VectorValue(x[1]^2,t)
∇f2(x,t) = ∇(y->f2(y,t))(x)
∇f2(xv,tv)
# @santiagobadia : Is there any way to make this transparent to the user
# I guess not unless we create a type for these analytical (space-only or
# space-time via a trait) functions
# Probably a try-catch?

# 2nd time derivative
f(x,t) = t^2
dtf = (x,t) -> ForwardDiff.derivative(t->f(x,t),t)
dttf = (x,t) -> ForwardDiff.derivative(t->dtf(x,t),t)
@test dttf(xv,tv) ≈ ∂tt(f)(xv,tv) ≈ ∂tt(f)(xv)(tv) ≈ ∂tt(f)(tv)(xv)
@test ∂tt(f)(xv,tv) ≈ 2.0

f2(x,t) = x[1]*t^2
dtf2 = (x,t) -> ForwardDiff.derivative(t->f2(x,t),t)
dttf2 = (x,t) -> ForwardDiff.derivative(t->dtf2(x,t),t)
@test dttf2(xv,tv) ≈ ∂tt(f2)(xv,tv) ≈ ∂tt(f2)(xv)(tv) ≈ ∂tt(f2)(tv)(xv)
@test ∂tt(f2)(xv,tv) ≈ 2.0*xv[1]

end #module
