# Trialize basis

function trialize(a::AbstractVector{<:Field})
  TrialFieldBasis(a)
end

struct TrialFieldBasis{B,T} <: AbstractArray{T,2}
  basis::B
  function TrialFieldBasis(basis)
    @assert ndims(basis) == 1
    T = eltype(basis)
    new{typeof(basis),T}(basis)
  end
end

Base.size(a::TrialFieldBasis) = (1,length(a.basis))
Base.axes(a::TrialFieldBasis) = (Base.OneTo(1),axes(a.basis,1)) # TODO OneTo
Base.getindex(a,TrialFieldBasis,i::Integer,j::Integer) = a.basis[j]
Base.getindex(a,TrialFieldBasis,i::Integer) = a.basis[i]
Base.IndexStyle(::TrialFieldBasis{A}) where A = IndexStyle(A)

function gradient(a::TrialFieldBasis)
  trialize(gradient(a))
end

function linear_combination(a::TrialFieldBasis,b::AbstractVector{<:Number})
  linear_combination(a.basis,b)
end

function return_cache(a::TrialFieldBasis,x)
  return_cache(a.basis,x)
end

function evaluate!(cache,a::TrialFieldBasis,x)
  M = evaluate!(cache,a.basis,x)
  trialize_on_values(M)
end

function trialize_on_values(M::AbstractVector)
  transpose(M)
end

function trialize_on_values(M::AbstractMatrix)
  TrializedMatrix(M)
end

struct TrializedMatrix{A,T} <: AbstractArray{T,3}
  matrix::A
  @inline function TrializedMatrix(matrix::AbstractMatrix{T}) where T
    A = typeof(matrix)
    new{A,T}(matrix)
  end
end

Base.size(a::TrializedMatrix) = (size(a.matrix,1),1,size(a.matrix,2))
Base.axes(a::TrializedMatrix) = (axes(a.matrix,1),Base.OneTo(1),axes(a.matrix,2)) # TODO
Base.IndexStyle(::Type{<:TrializedMatrix{A}}) where A = IndexStyle(A)
@inline Base.getindex(a::TrializedMatrix,i::Integer,j::Integer,k::Integer) = a.matrix[i,k]
@inline Base.getindex(a::TrializedMatrix,i::Integer) = a.matrix[i]
@inline Base.setindex!(a::TrializedMatrix,v,i::Integer,j::Integer,k::Integer) = (a.matrix[i,k] = v)
@inline Base.setindex!(a::TrializedMatrix,v,i::Integer) = (a.matrix[i] = v)
