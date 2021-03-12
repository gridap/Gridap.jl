function sparse_from_coo(::Type{<:SymSparseMatrixCSR{Bi}}, I,J,V,m,n) where Bi
  symsparsecsr(Val(Bi),M,I,J,V,m,n)
end

@inline function is_entry_stored(::Type{<:SymSparseMatrixCSR},i,j)
  (i<=j)
end

function finalize_coo!(T::Type{<:SymSparseMatrixCSR},I,J,V,m,n) 
  m == n || throw(DimensionMismatch("matrix is not square: dimensions are ($m, $n)"))
  touched = zeros(Bool,m)
  for k in 1:length(I)
    Ik = I[k]
    Jk = J[k]
    if Ik == Jk
      touched[Ik] = true
    end
  end
  for k in 1:m
    if ! touched[k]
      push!(I,k)
      push!(J,k)
      push!(V,zero(eltype(V)))
    end
  end
end

function nzindex(A::SparseMatrixCSR,i,j)
  if i<j
    nzindex(A.uppertrian,i,j)
  else
    nzindex(A.uppertrian,j,i)
  end
end

