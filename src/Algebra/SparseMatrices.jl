
# Extended Sparse matrix interface

"""
    sparse_from_coo(::Type{T} where T,I,J,V,m,n)
"""
function sparse_from_coo(::Type{T} where T,I,J,V,m,n)
  @abstractmethod
end

"""
    fill_entries!(A::AbstractSparseMatrix,v)

Fills the non-zero entries in the sparse matrix A with a given value v.
"""
function fill_entries!(A::AbstractSparseMatrix,v)
  @abstractmethod
end

"""
    push_coo!(::Type{T} where T, I,J,V,ik,jk,vk)

Inserts entries in COO vectors for further building a sparse matrix of type T.
"""
function push_coo!(::Type{T} where T,I,J,V,ik,jk,vk)
  @abstractmethod
end

"""
    is_entry_stored(::Type{T} where T,i,j) -> Bool

Tells if the entry with coordinates `[i,j]` will be stored when calling function push_coo!
"""
function is_entry_stored(::Type{T} where T,i,j)
  @abstractmethod
end


"""
    finalize_coo!(::Type{T} where T,I,J,V,m,n)

Check and insert diagonal entries in COO vectors if needed.
"""
function finalize_coo!(::Type{T} where T,I,J,V,m,n)
  @abstractmethod
end

"""
"""
function allocate_coo_vectors(::Type{M},n::Integer) where M
  return (zeros(Int,n), zeros(Int,n), zeros(Float64,n))
end

function allocate_coo_vectors(::Type{M},n::Integer) where {Tv,Ti,M<:AbstractSparseMatrix{Tv,Ti}}
  return (zeros(Ti,n), zeros(Ti,n), zeros(Tv,n))
end

