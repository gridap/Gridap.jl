
# SparseMatrices implementation contains:
# 
#   - Data types:
#     + `SparseMatrixCSR`:Compressed Sparse Row (CSR) sparse matrix implementation with Bi-based indexing.
#     + `SymSparseMatrixCSR`: Symmetric Compressed Sparse Row sparse matrix implementation with Bi-based indexing.
# 
#   - Procedures:
#     + `push_coo!`: Helper function to build COO arrays for further building a SparseMatrix
#     + `finalize_coo!`: Finalization of COO arrays building.
#     + `sparse_from_coo`: Return a SparseMatrix from COO data given.
#     + `add_entry!`: Add an entry given its position and the operation to perform.
#     + `fill_entries!`: Fills an existing sparse matrix with a given value.

# Extended Sparse matrix interface

"""
    sparse_from_coo(::Type{T} where T,args...)

`args...` are the same as for function `sparse`
"""
function sparse_from_coo(::Type{T} where T,args...)
  @abstractmethod
end

"""
    add_entry!(A,v::Number,i::Integer,j::Integer,combine::Function=+)

Add an entry given its position and the operation to perform.
"""
function add_entry!(A,v::Number,i::Integer,j::Integer,combine::Function=+)
  @abstractmethod
end

"""
    fill_entry!(A,v::Number)

Fills an existing sparse matrix A with a given value v.
"""
function fill_entries!(A,v::Number)
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
    is_entry_stored(::Type{T} where T,i::Integer,j::Integer) -> Bool

Tells if the entry with coordinates `[i,j]` will be stored when calling function push_coo!
"""
function is_entry_stored(::Type{T} where T,i::Integer,j::Integer)
  @abstractmethod
end


"""
    finalize_coo!(::Type{T} where T,I,J,V,m,n) 

Check and insert diagonal entries in COO vectors if needed.
"""
function finalize_coo!(::Type{T} where T,I::Vector,J::Vector,V::Vector,m::Integer,n::Integer)
  @abstractmethod
end

"""
"""
function create_coo_vectors(::Type{M}) where {M}
  return (Int[], Int[], Float64[])
end

function create_coo_vectors(::Type{M}) where {Tv,Ti,M<:AbstractSparseMatrix{Tv,Ti}}
  return (Ti[], Ti[], Tv[])
end

"""
"""
function allocate_coo_vectors(::Type{M},n::Integer) where M
  return (zeros(Int,n), zeros(Int,n), zeros(Float64,n))
end

function allocate_coo_vectors(::Type{M},n::Integer) where {Tv,Ti,M<:AbstractSparseMatrix{Tv,Ti}}
  return (zeros(Ti,n), zeros(Ti,n), zeros(Tv,n))
end

"""
"""
function copy_entries!(a::AbstractMatrix,b::AbstractMatrix)
  if a !== b
    _copy!(a,b)
  end
end

# We define this, since its not in 1.0
function _copy!(a,b)
  @assert size(a) == size(b) "Array dimension mismatch when copying"
  for i in eachindex(a)
    a[i] = b[i]
  end
end

