
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

function nz_counter(
  builder::SparseMatrixBuilder{SparseMatrixCSR{Bi,Tv,Ti}},axes) where {Bi,Tv,Ti}

  builder_csc = SparseMatrixBuilder(SparseMatrixCSC{Tv,Ti},builder.approach)
  csc =  nz_counter(builder_csc,(axes[2],axes[1]))
  bi = Val{Bi}()
  NzCounterCSR(bi,csc)
end

struct NzCounterCSR{Bi,T}
  bi::Val{Bi}
  csc::T
end

LoopStyle(::Type{NzCounterCSR{Bi,T}}) where {Bi,T} = LoopStyle(T)

@inline function add_entry!(::typeof(+),a::NzCounterCSR,v,i,j)
  add_entry!(+,a.csc,v,j,i)
end

function nz_allocation(a::NzCounterCSR)
  csc = nz_allocation(a.csc)
  NzAllocationCSR(a.bi,csc)
end

struct NzAllocationCSR{Bi,T}
  bi::Val{Bi}
  csc::T
end

LoopStyle(::Type{NzAllocationCSR{Bi,T}}) where {Bi,T} = LoopStyle(T)

@inline function add_entry!(::typeof(+),a::NzAllocationCSR,v,i,j)
  add_entry!(+,a.csc,v,j,i)
end

function create_from_nz(a::NzAllocationCSR{Bi}) where Bi
  Atcsc = create_from_nz(a.csc)
  Acsc = transpose(Atcsc)
  Acsr = SparseMatrixCSR{Bi}(Acsc)
  Acsr
end

