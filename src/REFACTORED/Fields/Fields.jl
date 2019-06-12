module Fields

using Test
using Gridap
using Gridap.Helpers

import Gridap: return_size
export FieldLike
export Field
export Basis
export Geomap
export num_dofs
export gradient
export ∇
export test_fieldlike
export test_field
export test_basis

"""
Umbrella type for Field and Basis
"""
const FieldLike{D,T<:FieldValue,N} = Map{Point{D},1,T,N}

"""
Create the gradient of a `Field` or `Basis`
For efficiency reasons, different calls to this function should return the same object
"""
function gradient(this::FieldLike{D,T,N})::FieldLike{D,G,N} where {D,T,G,N}
  @abstractmethod
end

const ∇ = gradient

"""
Abstract field of rank `T` (e.g., scalar, vector, tensor) on a manifold of
dimension `D`
"""
const Field{D,T} = FieldLike{D,T,1}

return_size(this::Field,s::Tuple{Int}) = s

"""
Abstract basis for a space of fields of rank `T` (e.g., scalar, vector, tensor)
on a manifold of dimension `D`
"""
const Basis{D,T} = FieldLike{D,T,2}

num_dofs(::Basis)::Int = @abstractmethod

function return_size(this::Basis,s::Tuple{Int})
  np, = s
  nd = num_dofs(this)
  (nd,np)
end

"""
Abstract geometry map
"""
const Geomap{D,Z,X} = Field{D,Point{Z,X}}

# Testers

function test_fieldlike(
  m::FieldLike{D,T,N},
  x::AbstractVector{<:Point{D}},
  v::AbstractArray{T,N},
  g::AbstractArray{G,N}) where {D,T,N,G}
  test_map(m,x,v)
  mg = gradient(m)
  test_map(mg,x,g)
  mg2 = gradient(m)
  @test mg === mg2
end

function test_field(
  m::Field{D,T},
  x::AbstractVector{<:Point{D}},
  v::AbstractVector{T},
  g::AbstractVector{G}) where {D,T,G}
  test_fieldlike(m,x,v,g)
end

function test_basis(
  m::Basis{D,T},
  x::AbstractVector{<:Point{D}},
  v::AbstractMatrix{T},
  g::AbstractMatrix{G}) where {D,T,G}
  nd = num_dofs(m)
  @test nd == size(v,1)
  @test nd == size(g,1)
  test_fieldlike(m,x,v,g)
  mg = gradient(m)
  nd = num_dofs(mg)
  @test nd == size(g,1)
end

end # module
