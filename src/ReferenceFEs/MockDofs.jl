
struct MockDofBasis{P} <: Dof
  x::Vector{P}
end

function dof_cache(b::MockDofBasis,field)
  field_cache(field,b.x)
end

@inline function evaluate_dof!(cache,b::MockDofBasis,field)
  evaluate_field!(cache,field,b.x)
end
