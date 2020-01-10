"""
    function push_coo!(::Type{SparseMatrixCSC},I,J,V,ik,jk,vk) 

Inserts entries in COO vectors for further building a SparseMatrixCSC.
"""
function push_coo!(::Type{SparseMatrixCSC},I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number)
    (push!(I, ik), push!(J, jk), push!(V, vk))
end

"""
    function push_coo!(I,J,V,ik,jk,vk) 

Inserts entries in COO vectors for further building a SparseMatrixCSC.
"""
function push_coo!(I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number)
    push_coo!(SparseMatrixCSC,I,J,V,ik,jk,vk)
end

"""
    function finalize_coo!(::Type{SparseMatrixCSC},I,J,V,m,n) 

Check and insert diagonal entries in COO vectors if needed.
"""
function finalize_coo!(::Type{SparseMatrixCSC},I::Vector,J::Vector,V::Vector,m::Integer,n::Integer)
end

"""
    function finalize_coo!(I,J,V,m,n) 

Finalize COO arrays for building a SparseMatrixCSC.
"""
function finalize_coo!(I::Vector,J::Vector,V::Vector,m::Integer,n::Integer)
    finalize_coo!(SparseMatrixCSC,I,J,V,m,n)
end


"""
    function hasrowmajororder(::Type{SparseMatrixCSC})

Check if values are stored in row-major order.
Return false.
"""
hasrowmajororder(::Type{SparseMatrixCSC}) = false
hasrowmajororder(a::SparseMatrixCSC) = hasrowmajororder(SparseMatrixCSC)

"""
    function hascolmajororder(::Type{SparseMatrixCSC})

Check if values are stored in col-major order.
Return true.
"""
hascolmajororder(::Type{SparseMatrixCSC}) = true
hascolmajororder(a::SparseMatrixCSC) = hascolmajororder(SparseMatrixCSC)

"""
    function getptr(S::SparseMatrixCSC)

Return columns pointer.
"""
getptr(S::SparseMatrixCSC) = S.colptr

"""
    function getvals(S::SparseMatrixCSC)

Return row indices.
"""
getindices(S::SparseMatrixCSC) = rowvals(S)


