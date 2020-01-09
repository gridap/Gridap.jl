
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

