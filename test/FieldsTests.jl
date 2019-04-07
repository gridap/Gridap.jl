module FieldsTests

##
using Test
using Numa.Fields
using Numa.FieldValues

import Numa.Fields: gradient
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
using StaticArrays
using Numa.FieldValues
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

end #module FieldsTests
