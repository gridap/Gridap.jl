
struct MockField{T,D} <: Field
  v::T
  function MockField{D}(v::Number) where {T,D}
    new{typeof(v),D}(v)
  end
end

MockField(D::Int,v::Number) = MockField{D}(v)

mock_field(D::Int,v::Number) = MockField{D}(v)

function evaluate!(c,f::MockField,x::Point)
  f.v*x[1]
end

@inline function evaluate_gradient!(cache,f::MockField{T,D},x::Point{D,S}) where {T,D,S}
  vg = outer(zero(Point{D,S}),zero(T))
end
