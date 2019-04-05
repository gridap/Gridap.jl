##
using Test
using Numa
using Numa.FieldValues
using Numa.Fields
using Numa.Polytopes

##

# Extension of Field for Analytical Fields
# Quite complicated... use Lambda-functions?
# Think how to be efficient for an array of points...
# Do we want to create evaluatefield! and evaluatefield?, etc
struct MyField{D,T} <: Field{D,T} end
import Numa.Fields: evaluatefield!
import Numa.Fields: evaluatefieldgradient!
function evaluatefield!(this::MyField{D,T}, points::Vector{Point{D}}, v::Vector{T}) where {D,T}
	for (p,P) in enumerate(points)
		v[p] = sum(P)*one(T)
	end
end
function evaluatefieldgradient!(this::MyField{D,T}, points::Vector{Point{D}}, v::Vector{TG}) where {D,T,TG}
	for (p,P) in enumerate(points)
		v[p] = ones(TG)
	end
end
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
