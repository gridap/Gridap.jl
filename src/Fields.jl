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
	MTG = mutable(eltype(v))
	z = zero(MTG)
	@inbounds for p in 1:length(points)
		z = zero(z)
		this.gradf(points[p],z)
		v[p] = z
	end
end

end  # module Fields
