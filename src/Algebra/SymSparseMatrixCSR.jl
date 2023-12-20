
function sparse_from_coo(::Type{<:SymSparseMatrixCSR{Bi}}, I,J,V,m,n) where Bi
  symsparsecsr(Val(Bi),I,J,V,m,n)
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

function nz_index(A::SymSparseMatrixCSR,i,j)
  nz_index(A.uppertrian,i,j)
end

function add_entry!(combine::Function,A::SymSparseMatrixCSR,v::Number,i,j)
  if i<=j
    k = nz_index(A,i,j)
    nz = nonzeros(A)
    Aij = nz[k]
    nz[k] = combine(Aij,v)
  end
  A
end

@inline function push_coo!(::Type{<:SymSparseMatrixCSR},I,J,V,i,j,v)
  (i>j) && return nothing
  push!(I,i)
  push!(J,j)
  push!(V,v)
  nothing
end

