function linear_combination(a::FieldVector, b::AbstractVector)
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
