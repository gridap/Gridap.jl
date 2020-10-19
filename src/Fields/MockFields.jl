
struct MockField{T<:Number} <: Field
  v::T
end

#MockField(D::Int,v::Number) = MockField{D}(v)

#mock_field(D::Int,v::Number) = MockField{D}(v)

function evaluate!(c,f::MockField,x::Point)
  f.v
end

function evaluate_gradient!(cache,f::MockField,x::Point)
  zero(outer(x,f.v))
end

function evaluate_hessian!(cache,f::MockField,x::Point)
  zero(outer(x,outer(x,f.v)))
end
