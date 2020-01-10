
# SparseMatrices implementation contains:
# 
#   - Data types:
#     + `SparseMatrixCSR`:Compressed Sparse Row (CSR) sparse matrix implementation with Bi-based indexing.
#     + `SymSparseMatrixCSR`: Symmetric Compressed Sparse Row sparse matrix implementation with Bi-based indexing.
# 
#   - Procedures:
#     + `push_coo!`: Helper function to build COO arrays for further building a SparseMatrix
#     + `finalize_coo!`: Finalization of COO arrays building.
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


import Base: convert, size, getindex, show, count, *
import LinearAlgebra: mul!
import SparseArrays: nnz, getnzval, nonzeros, nzrange, findnz, rowvals

# Extended Sparse matrix interface

"""
    sparse_from_coo(::Type{T} where T,args...)

`args...` are the same as for function `sparse`
"""
function sparse_from_coo(::Type{T} where T,args...)
  @abstractmethod
end

"""
    push_coo!(::Type{T} where T, I,J,V,ik,jk,vk)

Inserts entries in COO vectors for further building a sparse matrix of type T.
"""
function push_coo!(::Type{T} where T,I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number)
  @abstractmethod
end


"""
    finalize_coo!(::Type{T} where T,I,J,V,m,n) 

Check and insert diagonal entries in COO vectors if needed.
"""
function finalize_coo!(::Type{T} where T,I::Vector,J::Vector,V::Vector,m::Integer,n::Integer)
  @abstractmethod
end

# Implementation for SparseMatrixCSC

function sparse_from_coo(::Type{<:SparseMatrixCSC} where T,args...)
  sparse(args...)
end

@inline function push_coo!(
   ::Type{<:SparseMatrixCSC},I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number)
   (push!(I, ik), push!(J, jk), push!(V, vk))
end

function finalize_coo!(
  ::Type{<:SparseMatrixCSC},I::Vector,J::Vector,V::Vector,m::Integer,n::Integer)
  nothing
end

include("SparseMatrixCSC.jl")
include("SparseMatrixCSR.jl")
include("SymSparseMatrixCSR.jl")


