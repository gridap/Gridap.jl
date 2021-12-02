
struct MockDofBasis{P} <: AbstractVector{PointValue}
  nodes::Vector{P}
end

function return_cache(b::MockDofBasis,field)
  return_cache(field,b.x)
end

function evaluate!(cache,b::MockDofBasis,field)
  evaluate!(cache,field,b.x)
end

Base.size(a::MockDofBasis) = (length(a.nodes),)
Base.axes(a::MockDofBasis) = (axes(a.nodes,1),)
# @santiagobadia : Not sure we want to create the moment dofs
Base.getindex(a::MockDofBasis,i::Integer) = PointValue(a.x[i])
Base.IndexStyle(::MockDofBasis) = IndexLinear()
