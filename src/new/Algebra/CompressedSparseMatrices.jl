
# CompressedSparseMatrices implementation contains:
# 
#   - Data types:
#     + `SparseMatrixCSR`:Compressed Sparse Row (CSR) sparse matrix implementation with Bi-based indexing.
#     + `SymSparseMatrixCSR`: Symmetric Compressed Sparse Row sparse matrix implementation with Bi-based indexing.
# 
#   - Procedures:
#     + `sparsecsr`: Analogous method to [`sparse`](@ref) function to build a SparseMatrixCSR.
#     + `symsparsecsr`: Analogous method to [`sparse`](@ref) function to build a SymSparseMatrixCSR.
#     + `colvals`: Analogous method to [`rowvals`](@ref) to return `colvals` array.
#     + `hasrowmajororder`: Return `true` if matrix values are ordered by row.
#     + `hascolmajororder`: Return `true` if matrix values are ordered by column.
#     + `getptr`: Return the pointer array of a SparseMatrix (`rowptr` or `colptr` depending on the SparseMatrix type)
#     + `getindices`: Return the indices array of a SparseMatrix (`rowval` or `colval` depending on the SparseMatrix type)
# 
#   - Overloaded procedures:
#     + `*`: SparseMatrix-Vector product.
#     + `mul!`: SparseMatrix-Vector product.
#     + `nnz`: Return the number of stored (filled) elements in a sparse array.
#     + `nonzeros`: Return `nzval` array.
#     + `nzrange`: Return the range of indices for a particular row or column (Depending on the SparseMatrix type)
#     + `findnz`: Return a tuple (I, J, V) where I and J are the row and column indices.
#     + `rowvals`: Return row indices or raises an error (Depending on the SparseMatrix type)
#     + `convert`: Type conversion
#     + `count`: Count the number of elements in nonzeros(S) for which predicate pred returns true. 

import Base: convert, size, getindex, show, count, *
import LinearAlgebra: mul!
import SparseArrays: nnz, nonzeros, nzrange, findnz, rowvals

"""
    hasrowmajororder(::Type{AbstractSparseMatrix})

Check if values are stored in row-major order.
Return false.
"""
function hasrowmajororder(::Type{AbstractSparseMatrix}) 
  @abstractmethod
end

function hasrowmajororder(S::AbstractSparseMatrix)
  hasrowmajororder(typeof(S))
end

"""
    hascolmajororder(::Type{AbstractSparseMatrix})

Check if values are stored in col-major order.
Return true.
"""
function hascolmajororder(::Type{AbstractSparseMatrix})
  @abstractmethod
end

function hascolmajororder(S::AbstractSparseMatrix)
  hascolmajororder(typeof(S))
end

"""
    getptr(S::AbstractSparseMatrix)

Return columns pointer.
"""
function getptr(S::AbstractSparseMatrix)
  @abstractmethod
end

"""
    getindices(S::AbstractSparseMatrix)

Return row indices.
"""
function getindices(S::AbstractSparseMatrix)
  @abstractmethod
end

"""
    rowvals(S::AbstractSparseMatrix)

Return row indices or raises an error (Depending on the SparseMatrix type)
"""
function rowvals(S::AbstractSparseMatrix)
  @abstractmethod
end

"""
    colvals(S::AbstractSparseMatrix)

Return columns indices or raises an error (Depending on the SparseMatrix type)
"""
function colvals(S::AbstractSparseMatrix)
  @abstractmethod
end

"""
    nnz(S::AbstractSparseMatrix)

Returns the number of stored (filled) elements in a sparse array.
"""
function nnz(S::AbstractSparseMatrix)
  @abstractmethod
end

"""
    nonzeros(S::AbstractSparseMatrix)

Return a vector of the structural nonzero values in sparse array S. 
This includes zeros that are explicitly stored in the sparse array. 
The returned vector points directly to the internal nonzero storage of S, 
and any modifications to the returned vector will mutate S as well.
"""
function nonzeros(S::AbstractSparseMatrix)
  @abstractmethod
end

"""
    count(pred, S::AbstractSparseMatrix) -> Integer

Count the number of elements in nonzeros(S) for which predicate pred returns true. 
"""
function count(pred, S::AbstractSparseMatrix)
  @abstractmethod
end

"""
    nzrange(S::AbstractSparseMatrix, index::Integer) where {Bi}

Return the range of indices to the structural nonzero values of a 
sparse matrix index (Row or column depending on the compression type). 
"""
function nzrange(S::AbstractSparseMatrix, index::Integer)
  @abstractmethod
end

"""
    findnz(S::AbstractSparseMatrix)

Return a tuple (I, J, V) where I and J are the row and column indices 
of the stored ("structurally non-zero") values in sparse matrix A, 
and V is a vector of the values.
"""
function findnz(S::AbstractSparseMatrix)
  @abstractmethod
end

"""
    convert(::Type{AbstractSparseMatrix}, x::AbstractSparseMatrix)

Convert x to a value of the first type given.
"""
function convert(::Type{AbstractSparseMatrix}, x::AbstractSparseMatrix)
  @abstractmethod
end

include("SparseMatrixCSC.jl")
include("SparseMatrixCSR.jl")
include("SymSparseMatrixCSR.jl")
