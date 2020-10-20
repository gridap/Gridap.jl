
struct MockField{T<:Number} <: Field
  v::T
end

#MockField(D::Int,v::Number) = MockField{D}(v)

#mock_field(D::Int,v::Number) = MockField{D}(v)

function evaluate!(c,f::MockField,x::Point)
  f.v
end

function evaluate!(cache,f::FieldGradient{1,<:MockField},x::Point)
  zero(outer(x,f.object.v))
end

function evaluate!(cache,f::FieldGradient{2,<:MockField},x::Point)
  zero(outer(x,outer(x,f.object.v)))
end
