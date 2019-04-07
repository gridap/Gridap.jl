##
using Test
using Numa
using Numa.FieldValues
using Numa.Fields
using Numa.Polytopes

##

# Extension of Field for Analytical Fields
# struct MyField{D,T} <: Field{D,T} end
# import Numa.Fields: evaluatefield!
# import Numa.Fields: evaluatefieldgradient!
# function evaluatefield!(this::MyField{D,T}, points::Vector{Point{D}}, v::Vector{T}) where {D,T}
# 	for (p,P) in enumerate(points)
# 		v[p] = sum(P)*one(T)
# 	end
# end
# function evaluatefieldgradient!(this::MyField{D,T}, points::Vector{Point{D}}, v::Vector{TG}) where {D,T,TG}
# 	for (p,P) in enumerate(points)
# 		v[p] = ones(TG)
# 	end
# end
##
using Numa.Polytopes: PointInt
D = 2
nparts1d = 2
nparts = nparts1d*ones(Int64,D)
order=2
orders=order*ones(Int64,D)
extrusion = PointInt{D}(ones(Int64,D))
polytope = Polytope(extrusion)
nodes = NodesArray(polytope,orders)
##
D=2
typeof(nodes.coordinates[1]) <: Point{2}
p = Point{D}(1.0,1.0)
v = Numa.Fields.evaluatefield(MyField{D,ScalarValue}(),[p,p])
v
gv = Numa.Fields.evaluatefieldgradient(MyField{D,ScalarValue}(),[p,p])
@test v[1] = v[2] == 2.0
@test gv[1][1] == gv[1][2] == gv[2][1] == gv[2][2] == 1.0
##

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
a = Vector{VectorValue{D}}(undef,(l,))
MT = FieldValues.mutable(eltype(a))
p = Vector{Point{D}}(undef,(l,))
p = [ones(Point{D}) for i in 1:l]
println("++++++MVector++++++")
@time foo_mvector_point(a,p)
@time foo_mvector_point(a,p)
anfield = AnalyticalField{D,VectorValue{D}}(foo_mvector_point!,foo_mvector_point!)
@time v = Numa.Fields.evaluatefield(anfield,p)
@time Numa.Fields.evaluatefield!(anfield,p,v)
@time Numa.Fields.evaluatefield!(anfield,p,v)
# @santiagobadia : evaluatefield! allocates as much as
# evaluatefield but it does not happen when calling foo_mvector_point
# Commenting this.funct(points[p],z) or in evaluatefields! no allocations
# If I declare the function in Fields.jl it works too... somehow, we must
# "fix" the function interface in Fields.jl

##
