module Fields

using Numa #@fverdugo to be eliminated
using Numa.Helpers
using Numa.FieldValues

export Field

"""
Abstract field of rank `T` (e.g., scalar, vector, tensor) on a manifold of
dimension `D`
"""
abstract type Field{D,T} end
# @santiagobadia :  We must check whether we are constraining the environment
# space D, e.g., the point dimension in Point{D}, with the dimension of the
# manifold, e.g., the field value in VectorValue{MD}.

"""
Evaluate the field on a set of points
"""
function evaluate(this::Field{D,T},
	points::Vector{Point{D}})::Vector{T} where {D,T}
	@abstractmethod end

function evaluategradient(this::Field{D,T},
	points::Vector{Point{D}})::Vector{T} where {D,T}
	@abstractmethod end

# struct AnalyticalField{D,T} where {D,T}
# 	f
# 	gf
# end







end  # module Fields
