module Fields

using Test
using Gridap
using Gridap.Helpers

export HasGradientStyle
export GradientYesStyle
export GradientNotStyle

export FieldLike
export Field
export Basis
export Geomap
export num_dofs
export gradient
export ∇
export test_fieldlike
export test_field
export test_field_without_gradient
export test_basis
import Base: length

"""
Trait used to determine if a `CellField` and `CellBasis` type has gradient.
This trait is used to precompute gradients and store them for efficiency
"""
abstract type HasGradientStyle end

struct GradientYesStyle <: HasGradientStyle end

struct GradientNotStyle <: HasGradientStyle end

"""
Umbrella type for Field and Basis
"""
const FieldLike{D,T<:FieldValue,N} = Map{Point{D},1,T,N}

function HasGradientStyle(::T) where T <:FieldLike
  HasGradientStyle(T)
end

"""
Create the gradient of a `Field` or `Basis`
"""
function gradient(this::FieldLike)
  @abstractmethod
end

const ∇ = gradient

"""
Abstract field of rank `T` (e.g., scalar, vector, tensor) on a manifold of
dimension `D`
"""
const Field{D,T<:FieldValue} = FieldLike{D,T,1}

HasGradientStyle(::Type{<:Field}) = GradientNotStyle()

"""
Abstract basis for a space of fields of rank `T` (e.g., scalar, vector, tensor)
on a manifold of dimension `D`.

A Basis is evaluated at an array of Points and returns a matrix of values.
The first dimension in the returned matrix corresponds to the dofs of the basis,
whereas the second dimension corresponds to the evaluation points.
"""
const Basis{D,T<:FieldValue} = FieldLike{D,T,2}

HasGradientStyle(::Type{<:Basis}) = GradientYesStyle()

@inline function num_dofs(b::Basis)
  n, = return_size(b,(1,))
  n
end

@inline length(b::Basis) = num_dofs(b)

"""
Abstract geometry map
"""
const Geomap = Field{D,Point{Z,X}} where {D,Z,X}

# Testers

function test_fieldlike(
  m::FieldLike{D,T,N},
  x::AbstractVector{<:Point{D}},
  v::AbstractArray{T,N},
  g::AbstractArray{G,N}) where {D,T,N,G}
  test_map(m,x,v)
  mg = gradient(m)
  test_map(mg,x,g)
end

function test_field(
  m::Field{D,T},
  x::AbstractVector{<:Point{D}},
  v::AbstractVector{T},
  g::AbstractVector{G}) where {D,T,G}
  test_fieldlike(m,x,v,g)
end

function test_field_without_gradient(
  m::Field{D,T},
  x::AbstractVector{<:Point{D}},
  v::AbstractVector{T}) where {D,T}
  test_map(m,x,v)
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
  nd = length(mg)
  @test nd == size(g,1)
end

end # module
