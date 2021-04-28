
# Inverse fields

struct InverseField{F} <: Field
  original::F
end

inverse_map(a::Field) = InverseField(a)
inverse_map(a::InverseField) = a.original

function return_cache(a::InverseField,x::Point)
  return return_cache(a.original,x), return_cache(∇(a.original),x)
end

function evaluate!(caches,a::InverseField,x::Point)
  P = typeof(x)
  D = length(x)
  cache,∇cache = caches
  # Initial guess. (Can this be improved?)
  # We cannot use an SVector here; nlsolve expects a mutable object.
  # MVector also does not work (nlsolve then encounters nans).
  y₀ = [zero(P)...]
  # Function and its derivative
  f!(F,y) = F .= SVector{D}(Tuple(evaluate!(cache,a.original,P(y)) - x))
  j!(J,y) = J .= SMatrix{D,D}(Tuple(evaluate!(∇cache,∇(a.original),P(y))))
  # Solve
  res = nlsolve(f!,j!,y₀)
  @assert converged(res)
  # Extract solution
  y = res.zero
  return P(y)
end

function return_cache(a::InverseField,xs::AbstractVector{<:Point})
  return return_cache(a,testitem(xs))
end
function evaluate!(cache,a::InverseField,xs::AbstractVector{<:Point})
  map(x->evaluate!(cache,a,x), xs)
end
