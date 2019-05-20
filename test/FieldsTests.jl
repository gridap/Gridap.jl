# module FieldsTests

##
using Test
using Gridap.Maps
using Gridap.FieldValues

import Gridap: evaluate, gradient
##
fun(x::Point{2}) = x[1]*x[2] + x[1]
typeof(fun)

gradfun(x::Point{2}) = VectorValue(x[2]+1.0,x[1])

gradgradfun(x::Point{2}) = TensorValue(0.0,1.0,1.0,0.0)

gradient(::typeof(fun)) = gradfun

gradient(::typeof(gradfun)) = gradgradfun

p1 = VectorValue(0.1,1.0)
p2 = VectorValue(1.1,2.0)
p3 = VectorValue(1.4,5.0)

p = [p1,p2,p3]

f = AnalyticalField(fun,2)
@test evaluate(f,p) == fun.(p)

∇f = gradient(f)
@test evaluate(∇f,p) == gradfun.(p)

∇∇f = gradient(∇f)
@test evaluate(∇∇f,p) == gradgradfun.(p)


##
gradf = gradient(f)
valf = evaluate(f,p)
valgf = evaluate(gradf,p)
@test isa(f,Map)
g = -f
@test evaluate(g,p) == -1*valf
gg = -gradf
@test evaluate(gg,p) == -1*valgf
g = f+f
@test evaluate(g,p) == 2*valf
gg = gradf+gradf
@test evaluate(gg,p) == 2*valgf
g = f-f
@test evaluate(g,p) == 0*valf
g = inner(f,f)
@test evaluate(g,p) == valf.*valf
g = f*f
@test evaluate(g,p) == valf.*valf
g = inner(gradf,gradf)
@test evaluate(g,p) == broadcast(inner,valgf,valgf)
##
# Now tests with basis
using Gridap.Polytopes
using Gridap.RefFEs
using Gridap.Polytopes: PointInt
D = 2
orders=[1,1]
extrusion = PointInt{D}(1,1)
polytope = Polytope(extrusion)
reffe = LagrangianRefFE{D,ScalarValue}(polytope, orders)
bas = reffe.shfbasis
gradb = gradient(reffe.shfbasis)
@test isa(bas,Basis)
@test isa(bas,Map)
@test isa(gradb,Basis)
#
vals = [1.0,2.0,3.0,4.0]
using Gridap.Maps: FieldFromExpand
fef = FieldFromExpand(bas,vals)
@test evaluate(fef,p) ≈ evaluate(bas,p)'*vals
#
fun(x::Float64) = 2*x
using Gridap.Maps: FieldFromCompose
ff = FieldFromCompose(fun,f)
@test evaluate(ff,p) == evaluate(f,p)*2
#
gfun(x::Point{2}) = 2*x
g = AnalyticalField(gfun,2)
using Gridap.Maps: Geomap
isa(g,Geomap)
using Gridap.Maps: FieldFromComposeExtended
fce = FieldFromComposeExtended(fun,g,f)
pp = evaluate(g,p)
a1 = evaluate(f,pp)
a2 = broadcast(fun,a1)
@test evaluate(fce,p) == a2
##


##
using StaticArrays
using Gridap.FieldValues
function foo_mvector_point(a::Vector{T},p::Vector{Point{D}}) where {D,T}
  MT = FieldValues.mutable(eltype(a))
  z = zero(MT)
  @inbounds for i in 1:length(a)
    z = zero(z)
    foo_mvector_point!(p[i],z)
    a[i] = z
  end
end
function foo_mvector_point!(p,z)
  z[1] = p[1]+p[2]*2
  z[2] = p[2]-p[1]
end
##
##
l = 100000
D=2
p = [ones(Point{D}) for i in 1:l]
fun(x::Point{D}) = VectorValue(x[2]+1.0,x[1])
anfield = AnalyticalField(fun,D)
@time v = evaluate(anfield,p)
##

# end #module FieldsTests
