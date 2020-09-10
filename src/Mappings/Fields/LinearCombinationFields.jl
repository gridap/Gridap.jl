
"""
lincomb(a::Field,b::AbstractVector)

Returns a field obtained by the "linear combination" of
the value of the field basis `a` and the coefficient vector `b`.  The value of the resulting field
evaluated at a vector of points `x` is defined as

ax = evaluate(a,x)
ax*b

On the other hand, the gradient of the resulting field is defined as

∇ax = evaluate(gradient(a),x)
∇ax*b

"""
function linear_combination(a::NewField, b::AbstractVector)
  LinearCombinationField(a, b)
end

struct LinearCombinationField{A,B} <: NewField
  basis::A
  coefs::B
  @inline function LinearCombinationField(basis, coefs)
    A = typeof(basis)
    B = typeof(coefs)
    new{A,B}(basis, coefs)
  end
end

function return_cache(f::LinearCombinationField, x::AbstractArray{<:Point})
  ca = return_cache(f.basis, x)
  a = evaluate!(ca, f.basis, x)
  b = f.coefs
  ck = return_cache(MatVecMapping(), a, b)
  (ca, ck)
end

@inline function evaluate!(cache, f::LinearCombinationField, x::AbstractArray{<:Point})
  ca, ck = cache
  a = evaluate!(ca, f.basis, x)
  b = f.coefs
  evaluate!(ck, MatVecMapping(), a, b)
end

function gradient(f::LinearCombinationField)
  g = gradient(f.basis)
  LinearCombinationField(g,f.coefs)
end

# @santiagobadia : Here I must specialize the function in AppliedArray
# for optimization with arrays

# struct LinComValued <: Mapping end

# @inline function mapping_cache(k::LinComValued, a, b)
# LinComField(a,b)
# end

# @inline function apply_mapping!(f, k::LinComValued, a, b)
#   f.basis = a
#   f.coefs = b
#   f
# end

# function apply_gradient(k::LinComValued, a, b)
#   g = field_array_gradient(a)
#   LinearCombinationField(g, b)
# end

# function mapping_evaluate(k::LinComValued, x, a, b)
#   ax = evaluate_field_array(a, x)
#   apply(MulMapping(), ax, b)
# end
