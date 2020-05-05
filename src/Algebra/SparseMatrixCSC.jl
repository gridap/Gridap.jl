
function sparse_from_coo(::Type{<:SparseMatrixCSC}, args...)
  sparse(args...)
end

@inline function push_coo!(
   ::Type{<:SparseMatrixCSC},I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number)

   push!(I, ik)
   push!(J, jk)
   push!(V, vk)
   nothing
end

@inline function is_entry_stored(::Type{<:SparseMatrixCSC},i::Integer,j::Integer)
  true
end

function push_coo!(I::Vector,J::Vector,V::Vector,ik::Integer,jk::Integer,vk::Number)
  push_coo!(SparseMatrixCSC{Float64,Int},I,J,V,ik,jk,vk)
end

function finalize_coo!(
  ::Type{<:SparseMatrixCSC},I::Vector,J::Vector,V::Vector,m::Integer,n::Integer)
  nothing
end

function finalize_coo!(I::Vector,J::Vector,V::Vector,m::Integer,n::Integer)
    finalize_coo!(SparseMatrixCSC,I,J,V,m,n)
end

function add_entry!(A::SparseMatrixCSC{Tv,Ti},v::Number,i::Integer,j::Integer,combine::Function=+) where {Tv,Ti<:Integer}
    cv = convert(Tv, v)
    ci = convert(Ti, i)
    cj = convert(Ti, j)
    m,n = size(A)
    if !((1 <= ci <= m) & (1 <= cj <= n))
        throw(BoundsError(A, (ci,cj)))
    end
    cols = getptr(A)
    indices = getindices(A)
    vals = nonzeros(A)
    coljfirstk = Int(cols[cj])
    coljlastk = Int(cols[cj+1] - 1)
    searchk = searchsortedfirst(indices, ci, coljfirstk, coljlastk, Base.Order.Forward)
    if searchk <= coljlastk && indices[searchk] == ci
        # Column j contains entry A[i,j]. Update
        vold = vals[searchk]
        vals[searchk] = combine(vold, cv)
    else
        throw(ArgumentError("Entry ($i,$j) does not exist in the sparsity pattern"))
    end
    return A
end

function copy_entries!(a::SparseMatrixCSC,b::SparseMatrixCSC)
  _copy_entries_sparse!(a,b)
end

function _copy_entries_sparse!(a,b)
  na = nonzeros(a)
  nb = nonzeros(b)
  if na !== nb
    copyto!(na,nb)
  end
end

hasrowmajororder(::Type{<:SparseMatrixCSC}) = false
hasrowmajororder(a::SparseMatrixCSC) = hasrowmajororder(SparseMatrixCSC)

hascolmajororder(::Type{<:SparseMatrixCSC}) = true
hascolmajororder(a::SparseMatrixCSC) = hascolmajororder(SparseMatrixCSC)

getptr(S::SparseMatrixCSC) = S.colptr

getindices(S::SparseMatrixCSC) = rowvals(S)


