
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
#     + `fill_entry!`: Fills an existing sparse matrix with a given value.

# Extended Sparse matrix interface

"""
    sparse_from_coo(::Type{T} where T,args...)

`args...` are the same as for function `sparse`
"""
function sparse_from_coo(::Type{T} where T,args...)
  @abstractmethod
end

"""
    add_entry!(A,v::Number,i::Integer,j::Integer,,combine::Function=+)

Add an entry given its position and the operation to perform.
"""
function add_entry!(A,v::Number,i::Integer,j::Integer,combine::Function=+)
  @abstractmethod
end

"""
    fill_entry!(A,v::Number)

Fills an existing sparse matrix A with a given value v.
"""
function fill_entry!(A,v::Number)
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


