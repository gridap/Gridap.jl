
# Inverse fields

struct InverseField{F} <: Field
  original::F
end

inverse_map(a::Field) = InverseField(a)
inverse_map(a::InverseField) = a.original

function return_cache(a::InverseField,x::Point)
  y₀ = [zero(x)...]             # initial guess and solution
  F₀ = [zero(x)...]             # error
  return return_cache(a.original,x), return_cache(∇(a.original),x), y₀, F₀
end

function evaluate!(caches,a::InverseField,x::Point)
  P = typeof(x)
  D = length(x)
  cache,∇cache,y₀,F₀ = caches
  # Initial guess
  # We cannot use an SVector here; nlsolve expects a mutable object.
  # We're keeping the vector in the cache to avoid allocations.
  y₀ .= Tuple(zero(x))
  # Function and its derivative
  f!(F,y) = F .= SVector{D}(Tuple(evaluate!(cache,a.original,P(y)) - x))
  j!(J,y) = J .= SMatrix{D,D}(Tuple(evaluate!(∇cache,∇(a.original),P(y))))
  df = OnceDifferentiable(f!,j!,y₀,F₀)
  # Solve
  res = nlsolve(df,y₀)
  @check converged(res) "InverseField evaluation did not converge"
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
