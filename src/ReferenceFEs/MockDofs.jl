
struct MockDofBasis{P} <: AbstractVector{PointValue}
  nodes::Vector{P}
end

function return_cache(b::MockDofBasis,field)
  return_cache(field,b.x)
end

@inline function evaluate!(cache,b::MockDofBasis,field)
  evaluate!(cache,field,b.x)
end

@inline Base.size(a::MockDofBasis) = (length(a.nodes),)
@inline Base.axes(a::MockDofBasis) = (axes(a.nodes,1),)
# @santiagobadia : Not sure we want to create the moment dofs
@inline Base.getindex(a::MockDofBasis,i::Integer) = PointValue(a.x[i])
@inline Base.IndexStyle(::MockDofBasis) = IndexLinear()
