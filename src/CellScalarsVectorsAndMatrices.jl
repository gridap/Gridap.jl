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

inner(test::CellFieldValues,trial::CellFieldValues) = CellScalarsFromInner(test,trial)
