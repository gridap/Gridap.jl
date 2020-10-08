
struct MockDofBasis{P} <: Dof
  x::Vector{P}
end

function dof_cache(b::MockDofBasis,field)
  return_cache(field,b.x)
end

@inline function evaluate_dof!(cache,b::MockDofBasis,field)
  evaluate!(cache,field,b.x)
end
