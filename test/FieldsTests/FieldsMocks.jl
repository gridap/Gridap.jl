module FieldsMocks

using Gridap
using Gridap.Helpers
using TensorValues
using StaticArrays

export MockField
export MockGeomap
export MockBasis

import Gridap: evaluate!
import Gridap: return_size
import Gridap: gradient

struct GradMockField{D,X} <: Field{D,VectorValue{D,X}} end

struct MockField{D,X} <: Field{D,X}
  g::GradMockField{D,X}
end

function MockField(D::Int,X::Type)
  g = GradMockField{D,X}()
  MockField{D,X}(g)
end

function evaluate!(
  this::MockField{D,X},
  points::AbstractVector{<:Point{D}},
  v::AbstractVector{X}) where {D,X}
  for i in eachindex(points)
    p = points[i]
    v[i] = p[1]
  end
end

gradient(f::MockField) = f.g

return_size(::MockField,s::NTuple{N,Int} where N) = s

return_size(::GradMockField,s::NTuple{N,Int} where N) = s

function evaluate!(
  this::GradMockField{D,X},
  points::AbstractVector{<:Point{D}},
  v::AbstractVector{VectorValue{D,X}}) where {D,X}
  z = zero(MVector{D,X})
  z[1] = one(X)
  for i in eachindex(points)
    v[i] = z
  end
end

gradient(f::GradMockField) = @notimplemented

struct GradMockBasis{D,X} <: Basis{D,VectorValue{D,X}} end

struct MockBasis{D,X} <: Basis{D,X}
  g::GradMockBasis{D,X}
end

function MockBasis(D::Int,X::Type)
  g = GradMockBasis{D,X}()
  MockBasis{D,X}(g)
end

function evaluate!(
  this::MockBasis{D,X},
  points::AbstractVector{<:Point{D}},
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
  points::AbstractVector{<:Point{D}},
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

function return_size(::MockBasis,s::NTuple{N,Int} where N)
  n, = s
  (3,n)
end

function return_size(::GradMockBasis,s::NTuple{N,Int} where N)
  n, = s
  (3,n)
end

gradient(f::MockBasis) = f.g

gradient(f::GradMockBasis) = @notimplemented

struct GradMockGeomap{D,X,L} <: Field{D,TensorValue{D,X,L}} end

struct MockGeomap{D,X} <: Field{D,Point{D,X}}
  g::GradMockGeomap{D,X}
end

function MockGeomap(D::Int,X::Type)
  g = GradMockGeomap{D,X,D*D}()
  MockGeomap{D,X}(g)
end

function evaluate!(
  this::MockGeomap{D,X},
  points::AbstractVector{<:Point{D}},
  v::AbstractVector{VectorValue{D,X}}) where {D,X}
  for i in eachindex(points)
    p = points[i]
    v[i] = 3*p
  end
end

gradient(f::MockGeomap) = f.g

return_size(::MockGeomap,s::NTuple{N,Int} where N) = s

return_size(::GradMockGeomap,s::NTuple{N,Int} where N) = s

function evaluate!(
  this::GradMockGeomap{D,X,L},
  points::AbstractVector{<:Point{D}},
  v::AbstractVector{TensorValue{D,X,L}}) where {D,X,L}
  z = 3*one(TensorValue{D,X,L})
  for i in eachindex(points)
    v[i] = z
  end
end

gradient(f::GradMockGeomap) = @notimplemented

end # module
