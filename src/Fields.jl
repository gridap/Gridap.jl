module Fields

using Numa #@fverdugo to be eliminated
using Numa.Helpers
using Numa.FieldValues

export Field
export AnalyticalField
export evaluatefield
export evaluategradient
export evaluatefield!
export evaluategradient!

"""
Abstract field of rank `T` (e.g., scalar, vector, tensor) on a manifold of
dimension `D`
"""
abstract type Field{D,T} end
# @santiagobadia :  We must check whether we are constraining the environment
# space D, e.g., the point dimension in Point{D}, with the dimension of the
# manifold, e.g., the field value in VectorValue{MD}.

# @fverdugo I would use `evaluate!` instead of `evaluatefield!` (julia dispatch will do the rest)
# Idem for related functions `evaluatefieldgradients!`, `evaluatefield`, etc.
# I also would define the abstract interface gradient(::Field)::Field = @abstractmethod in order to
# be consistent with other parts of the code, which in turn will allow to reduce a lot of code replication
# See Fields2.jl as an example of the proposed interface

#@fverdugo v::AbstractVector{T}
"""
Evaluate the field on a set of points
"""
function evaluatefield!(this::Field{D,T},
	points::AbstractVector{Point{D}}, v::Vector{T}) where {D,T}
	@abstractmethod
end

function evaluatefieldgradient!(this::Field{D,T},
	points::AbstractVector{Point{D}}, v::Vector{TG}) where {D,T,TG}
	@abstractmethod
end


"""
Same as evaluate! but allocates output
"""
function evaluatefield(this::Field{D,T},
	points::AbstractVector{Point{D}}) where {D,T}
  vals = Vector{T}(undef,length(points))
  evaluatefield!(this, points, vals)
  vals
end

function evaluatefieldgradient(this::Field{D,T},
	points::AbstractVector{Point{D}}) where {D,T}
	TG = outer(Point{D},T)
  vals = Vector{TG}(undef,length(points))
  evaluatefieldgradient!(this, points, vals)
  vals
end

# @fverdugo this type is type unstable, and it is limited to 2nd other derivatives.
# See file Fields2.jl and the api in FieldsTests2.jl for an alternative implementation
# which fixes both issues.
"""
Field generated from user-provided lambda-functions for field value and its
gradient
"""
struct AnalyticalField{D,T} <: Field{D,T}
	funct
	gradf
	# @santiagobadia : Even this value is exported with analyticfield and creates
	# conflicts in other tests that define a f function
end
# @santiagobadia: How can we enforce f : Point{D} -> T and
# g : Point{D} -> TG
function evaluatefield!(this::AnalyticalField{D,T}, points::Vector{Point{D}}, v::Vector{T}) where {D,T}
	MT = mutable(eltype(v))
	z = zero(MT)
	@inbounds for p in 1:length(points)
		z = zero(z)
		this.funct(points[p],z)
		v[p] = z
	end
end

function evaluatefieldgradient!(this::AnalyticalField{D,T}, points::Vector{Point{D}}, v::Vector{TG}) where {D,T,TG}
	for (p,P) in enumerate(points)
		v[p] = this.gradf(P)
	end
end
# struct AnalyticalField{D,T} where {D,T}
# 	f
# 	gf
# end

# struct MyField{D,T} <: Field{D,T} end
#
# function evaluatefield(this::MyField{D,T}, point::Point{D})
# 	return sum(point)*one(T)
# end
#
# function evaluategradient(this::MyField{D,T}, point::Point{D})
#   TG = outer(Point{D},T)
# 	return one(TG)
# end







end  # module Fields
