
function sparse_from_coo(::Type{<:SparseMatrixCSC} where T,args...)
  sparse(args...)
end

@inline function push_coo!(
   ::Type{<:SparseMatrixCSC},I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number)
   (push!(I, ik), push!(J, jk), push!(V, vk))
end

function push_coo!(I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number)
    push_coo!(SparseMatrixCSC,I,J,V,ik,jk,vk)
end

function finalize_coo!(
  ::Type{<:SparseMatrixCSC},I::Vector,J::Vector,V::Vector,m::Integer,n::Integer)
  nothing
end

function finalize_coo!(I::Vector,J::Vector,V::Vector,m::Integer,n::Integer)
    finalize_coo!(SparseMatrixCSC,I,J,V,m,n)
end

hasrowmajororder(::Type{<:SparseMatrixCSC}) = false
hasrowmajororder(a::SparseMatrixCSC) = hasrowmajororder(SparseMatrixCSC)

hascolmajororder(::Type{<:SparseMatrixCSC}) = true
hascolmajororder(a::SparseMatrixCSC) = hascolmajororder(SparseMatrixCSC)

getptr(S::SparseMatrixCSC) = S.colptr

getindices(S::SparseMatrixCSC) = rowvals(S)


