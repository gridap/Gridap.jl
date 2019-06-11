module Fields

using Gridap
using Gridap.Helpers

import Gridap: evaluate!
import Gridap: return_size
export FieldLike
export Field
export Basis
export Geomap
export gradient
export ∇
export test_fieldlike
export test_field
export test_basis

"""
Umbrella type for Field and Basis
"""
const FieldLike{D,T<:FieldValue,N,X} = Map{Point{D,X},1,T,N}

"""
Create the gradient of a `Field` or `Basis`
"""
function gradient(this::FieldLike{D,T,N})::FieldLike{D,G,N} where {D,T,G,N}
  @abstractmethod
end

const ∇ = gradient

"""
Abstract field of rank `T` (e.g., scalar, vector, tensor) on a manifold of
dimension `D`
"""
const Field{D,T,X} = FieldLike{D,T,1,X}

"""
Abstract basis for a space of fields of rank `T` (e.g., scalar, vector, tensor)
on a manifold of dimension `D`
"""
const Basis{D,T,X} = FieldLike{D,T,2,X}

"""
Abstract geometry map
"""
const Geomap{D,Z,X} = Field{D,Point{Z,X},X}

# Testers

function test_fieldlike(
  m::FieldLike{D,T,N,X},
  x::AbstractVector{Point{D,X}},
  v::AbstractArray{T,N},
  g::AbstractArray{G,N}) where {D,T,N,G,X}
  test_map(m,x,v)
  mg = gradient(m)
  test_map(mg,x,g)
end

function test_field(
  m::Field{D,T,X},
  x::AbstractVector{Point{D,X}},
  v::AbstractVector{T},
  g::AbstractVector{G}) where {D,T,G,X}
  test_fieldlike(m,x,v,g)
end

function test_basis(
  m::Basis{D,T,X},
  x::AbstractVector{Point{D,X}},
  v::AbstractMatrix{T},
  g::AbstractMatrix{G}) where {D,T,G,X}
  test_fieldlike(m,x,v,g)
end

end # module
