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

function inner(test::CellFieldValues{T},trial::CellFieldValues{T}) where T

  function computesize(test,trial)
    @assert test == trial
    test
  end

  function computevals!(test,trial,values)
    @assert length(test) <= length(values)
    @assert length(trial) <= length(values)
    npoin = length(trial)
    @inbounds for i = 1:npoin
      values[i] = inner(test[i],trial[i])
    end
    @view values[1:npoin]
  end

  CellArrayFromBinaryOp{T,1,T,1,Float64,1}(test,trial,computevals!,computesize)

end

function inner(test::CellBasisValues{T},trial::CellFieldValues{T}) where T

  function computesize(test,trial)
    @assert test[2] == trial[1]
    test
  end

  function computevals!(test,trial,values)
    @assert size(test,2) <= length(values)
    @assert size(test,2) == length(trial)
    ndofs = size(test,1)
    npoin = size(test,2)
    @inbounds for i = 1:npoin
      @inbounds for j = 1:ndofs
        values[j,i] = inner(test[j,i],trial[i])
      end
    end
    @view values[1:ndofs,1:npoin]
  end

  CellArrayFromBinaryOp{T,2,T,1,Float64,2}(test,trial,computevals!,computesize)

end

function inner(test::CellBasisValues{T},trial::CellBasisValues{T}) where T

  function computesize(test,trial)
    @assert test[2] == trial[2]
    (test[1],trial[1],test[2])
  end

  function computevals!(test,trial,values)
    @assert size(test,2) <= size(values,3)
    @assert size(test,1) <= size(values,1)
    @assert size(trial,2) <= size(values,3)
    @assert size(trial,1) <= size(values,2)
    npoints = size(test,2)
    ndofstest = size(test,1)
    ndofstrial = size(trial,1)
    @inbounds for i = 1:npoints
      @inbounds for j = 1:ndofstrial
        @inbounds for k = 1:ndofstest
          values[k,j,i] = inner(test[k,i],trial[j,i])
        end
      end
    end
    @view values[1:ndofstest,1:ndofstrial,1:npoints]
  end

  CellArrayFromBinaryOp{T,2,T,2,Float64,3}(test,trial,computevals!,computesize)

end

function inner(test::CellField{D,T},trial::CellField{D,T}) where {D,T}
  EvaluableCellArrayFromBinaryOp{D,T,1,1,1}(test,trial,inner)
end

function inner(test::CellBasis{D,T},trial::CellField{D,T}) where {D,T}
  EvaluableCellArrayFromBinaryOp{D,T,2,2,1}(test,trial,inner)
end

function inner(test::CellBasis{D,T},trial::CellBasis{D,T}) where {D,T}
  EvaluableCellArrayFromBinaryOp{D,T,3,2,2}(test,trial,inner)
end

