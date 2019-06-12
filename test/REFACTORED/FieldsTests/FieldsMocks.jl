module FieldsMocks

using Gridap
using Gridap.Helpers
using TensorValues
using StaticArrays

export MockField
export MockBasis

import Gridap: evaluate!
import Gridap: return_size
import Gridap: num_dofs
import Gridap: gradient

struct GradMockField{D,X} <: Field{D,VectorValue{D,X},X} end

struct MockField{D,X} <: Field{D,X,X}
  g::GradMockField{D,X}
end

function MockField(D::Int,X::Type)
  g = GradMockField{D,X}()
  MockField{D,X}(g)
end

function evaluate!(
  this::MockField{D,X},
  points::AbstractVector{Point{D,X}},
  v::AbstractVector{X}) where {D,X}
  for i in eachindex(points)
    p = points[i]
    v[i] = p[1]
  end
end

gradient(f::MockField) = f.g

function evaluate!(
  this::GradMockField{D,X},
  points::AbstractVector{Point{D,X}},
  v::AbstractVector{VectorValue{D,X}}) where {D,X}
  z = zero(MVector{D,X})
  z[1] = one(X)
  for i in eachindex(points)
    v[i] = z
  end
end

gradient(f::GradMockField) = @notimplemented

struct GradMockBasis{D,X} <: Basis{D,VectorValue{D,X},X} end

struct MockBasis{D,X} <: Basis{D,X,X}
  g::GradMockBasis{D,X}
end

function MockBasis(D::Int,X::Type)
  g = GradMockBasis{D,X}()
  MockBasis{D,X}(g)
end

function evaluate!(
  this::MockBasis{D,X},
  points::AbstractVector{Point{D,X}},
  v::AbstractMatrix{X}) where {D,X}
  for j in eachindex(points)
    p = points[j]
    for i in 1:3
      v[i,j] = i*p[1]
    end
  end
end

function evaluate!(
  this::GradMockBasis{D,X},
  points::AbstractVector{Point{D,X}},
  v::AbstractMatrix{VectorValue{D,X}}) where {D,X}
  z = zero(MVector{D,X})
  z[1] = one(X)
  for j in eachindex(points)
    p = points[j]
    for i in 1:3
      v[i,j] = i*z
    end
  end
end

num_dofs(::MockBasis) = 3

num_dofs(::GradMockBasis) = 3

gradient(f::MockBasis) = f.g

gradient(f::GradMockBasis) = @notimplemented

end # module
