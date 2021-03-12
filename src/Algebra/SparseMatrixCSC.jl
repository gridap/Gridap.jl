
function sparse_from_coo(::Type{<:SparseMatrixCSC},I,J,V,m,n)
  sparse(I,J,V,m,n)
end

@inline function is_entry_stored(::Type{<:SparseMatrixCSC},i,j)
  true
end

function finalize_coo!(::Type{<:SparseMatrixCSC},I,J,V,m,n)
  nothing
end

function nz_index(A::SparseArrays.AbstractSparseMatrixCSC,i0,i1)
    if !(1 <= i0 <= size(A, 1) && 1 <= i1 <= size(A, 2)); throw(BoundsError()); end
    ptrs = SparseArrays.getcolptr(A)
    r1 = Int(ptrs[i1])
    r2 = Int(ptrs[i1+1]-1)
    (r1 > r2) && return -1
    r1 = searchsortedfirst(rowvals(A), i0, r1, r2, Base.Order.Forward)
    ((r1 > r2) || (rowvals(A)[r1] != i0)) ? -1 : r1
end

@inline function push_coo!(::Type{<:SparseMatrixCSC},I,J,V,i,j,v)
 push!(I,i)
 push!(J,j)
 push!(V,v)
 nothing
end
