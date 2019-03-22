export inner
"""
Abstract type that represents a scalar of value T
associated with a collection of points in each cell
"""
const CellScalars{T} = CellArray{T,1} where T

"""
Abstract type that represents a vector of value T
associated with a collection of points in each cell
(typically the cell rhs vector at the quadrature points)
"""
const CellVectors{T} = CellArray{T,2} where T

"""
Abstract type that represents a matrix of value T
associated with a collection of points in each cell
(typically the cell matrix at the quadrature points)
"""
const CellMatrices{T} = CellArray{T,3} where T


struct CellScalarsFromInner{T} <: CellScalars{T}
  test::CellFieldValues{T}
  trial::CellFieldValues{T}
end

function Base.length(self::CellScalarsFromInner)
  @assert length(self.test) == length(self.trial)
  length(self.test)
end

function maxsize(self::CellScalarsFromInner)
  @assert maxsize(self.test) == maxsize(self.trial)
  maxsize(self.test)
end

function Base.iterate(self::CellScalarsFromInner{T}) where {T}
  maxqpoints = maxsize(self.test,1)
  values = Array{T,1}(undef,(maxqpoints,))
  testnext = iterate(self.test)
  trialnext = iterate(self.trial)
  state = (values,testnext,trialnext)
  iterate(self,state)
end

function Base.iterate(self::CellScalarsFromInner,state)
  (values,testnext,trialnext) = state
  if testnext == nothing || trialnext == nothing
    nothing
  else
    (testvals,teststate) = testnext
    (trialvals,trialstate) = trialnext
    inner_scalar_kernel!(testvals,trialvals,values)
    npoints = size(testvals,1)
    vals = @view values[1:npoints]
    testnext = iterate(self.test,teststate)
    trialnext = iterate(self.trial,trialstate)
    state = (values,testnext,trialnext)
    (vals, state)
  end
end

function inner_scalar_kernel!(testvals,trialvals,values)
  @assert length(testvals) <= length(values)
  @assert length(trialvals) <= length(values)
  npoin = length(trialvals)
  @inbounds for i = 1:npoin
    values[i] = inner(testvals[i],trialvals[i])
  end
end

#TODO this involves a lot of code repetition
# Surely there is a better way to do it in julia
# We want DRY (don't repeat yourself) code

struct CellVectorsFromInner{T} <: CellVectors{T}
  test::CellBasisValues{T}
  trial::CellFieldValues{T}
end

function Base.length(self::CellVectorsFromInner)
  @assert length(self.test) == length(self.trial)
  length(self.test)
end

function maxsize(self::CellVectorsFromInner)
  @assert maxsize(self.test,2) == maxsize(self.trial,1)
  maxsize(self.test)
end

function Base.iterate(self::CellVectorsFromInner{T}) where {T}
  maxdofs = maxsize(self.test,1)
  maxqpoints = maxsize(self.test,2)
  values = Array{T,2}(undef,(maxdofs,maxqpoints))
  testnext = iterate(self.test)
  trialnext = iterate(self.trial)
  state = (values,testnext,trialnext)
  iterate(self,state)
end

function Base.iterate(self::CellVectorsFromInner,state)
  (values,testnext,trialnext) = state
  if testnext == nothing || trialnext == nothing
    nothing
  else
    (testvals,teststate) = testnext
    (trialvals,trialstate) = trialnext
    inner_vector_kernel!(testvals,trialvals,values)
    ndofs = size(testvals,1)
    npoints = size(testvals,2)
    vals = @view values[1:ndofs,1:npoints]
    testnext = iterate(self.test,teststate)
    trialnext = iterate(self.trial,trialstate)
    state = (values,testnext,trialnext)
    (vals, state)
  end
end

function inner_vector_kernel!(testvals,trialvals,values)
  @assert size(testvals,2) <= length(values)
  @assert length(trialvals) <= length(values)
  ndofs = size(testvals,1)
  npoin = size(testvals,2)
  @inbounds for i = 1:npoin
    @inbounds for j = 1:ndofs
      values[j,i] = inner(testvals[j,i],trialvals[i])
    end
  end
end

inner(test::CellFieldValues,trial::CellFieldValues) = CellScalarsFromInner(test,trial)

inner(test::CellBasisValues,trial::CellFieldValues) = CellVectorsFromInner(test,trial)
