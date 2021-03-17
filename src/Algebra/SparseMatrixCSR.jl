
function sparse_from_coo(::Type{<:SparseMatrixCSR{Bi}}, I,J,V,m,n) where Bi
  sparsecsr(Val(Bi),I,J,V,m,n)
end

@inline function is_entry_stored(::Type{<:SparseMatrixCSR},i,j)
  true
end

function finalize_coo!(::Type{<:SparseMatrixCSR},I,J,V,m,n)
  nothing
end

function nz_index(A::SparseMatrixCSR,i0,i1)
  if !(1 <= i0 <= size(A, 1) && 1 <= i1 <= size(A, 2)); throw(BoundsError()); end
  o = getoffset(A)
  Bi = getBi(A)
  r1 = Int(A.rowptr[i0]+o)
  r2 = Int(A.rowptr[i0+1]-Bi)
  (r1 > r2) && return -1
  i1o = i1-o
  k = searchsortedfirst(colvals(A), i1o, r1, r2, Base.Order.Forward)
  ((k > r2) || (colvals(A)[k] != i1o)) ? -1 : k
end

@inline function push_coo!(::Type{<:SparseMatrixCSR},I,J,V,i,j,v)
 push!(I,i)
 push!(J,j)
 push!(V,v)
 nothing
end
