
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

#@inline function add_entries!(
#  combine::Function,
#  A::SparseMatrixCSC,
#  vs::AbstractMatrix{<:Number},
#  is,js)
#
#  if issorted(is)
#    nz = A.nzval
#    ptrs = A.colptr
#    rows = A.rowval
#    for (lj,j) in enumerate(js)
#      if j>0
#        pini = ptrs[j]
#        pend = ptrs[j+1]-1
#        li = 1
#        for p in pini:pend
#          _i = rows[p]
#          if _i == is[li]
#            vij = vs[li,lj]
#            Aij = nz[p]
#            nz[p] = combine(Aij,vij)
#            li += 1
#          end
#        end
#      end
#    end
#  else
#    for (lj,j) in enumerate(js)
#      if j>0
#        for (li,i) in enumerate(is)
#          if i>0
#            vij = vs[li,lj]
#            add_entry!(combine,A,vij,i,j)
#          end
#        end
#      end
#    end
#  end
#  A
#end

struct CounterCSRR{Tv,Ti}
  tv::Type{Tv}
  nrows::Int
  ncols::Int
  rowptrs::Vector{Ti}
end

LoopStyle(::Type{<:CounterCSRR}) = Loop()

@inline function add_entry!(::typeof(+),a::CounterCSRR{Tv,Ti},v,i,j) where {Tv,Ti}
  a.rowptrs[i+1] += Ti(1)
  nothing
end

struct CSRR{Tv,Ti}
  nrows::Int
  ncols::Int
  rowptrs::Vector{Ti}
  colvals::Vector{Ti}
  nzvals::Vector{Tv}
end

LoopStyle(::Type{<:CSRR}) = Loop()

@inline function add_entry!(::typeof(+),a::CSRR{Tv,Ti},v::Nothing,i,j) where {Tv,Ti}
  p = a.rowptrs[i]
  a.colvals[p] = j
  a.rowptrs[i] = p+Ti(1)
  nothing
end

@inline function add_entry!(::typeof(+),a::CSRR{Tv,Ti},v,i,j) where {Tv,Ti}
  p = a.rowptrs[i]
  a.colvals[p] = j
  a.nzvals[p] = v
  a.rowptrs[i] = p+Ti(1)
  nothing
end

function nz_counter(::Type{SparseMatrixCSC{Tv,Ti}},axes) where {Tv,Ti}
  nrows = length(axes[1])
  ncols = length(axes[2])
  rowptrs = zeros(Ti,nrows+1)
  CounterCSRR(Tv,nrows,ncols,rowptrs)
end

function nz_allocation(a::CounterCSRR{Tv,Ti}) where {Tv,Ti}
  rowptrs = a.rowptrs
  length_to_ptrs!(rowptrs)
  ndata = rowptrs[end]-1
  colvals = zeros(Ti,ndata)
  nzvals = zeros(Tv,ndata)
  CSRR(a.nrows,a.ncols,rowptrs,colvals,nzvals)
end

function create_from_nz(a::CSRR{Tv,Ti}) where {Tv,Ti}
  rewind_ptrs!(a.rowptrs)
  A = _csrr_to_csc!(a)
  A
end

function _csrr_to_csc!(csrr::CSRR{Tv,Ti}) where {Tv,Ti}
  nrows = csrr.nrows
  ncols = csrr.ncols
  rowptrs = csrr.rowptrs
  colvals = csrr.colvals
  nzvalscsr = csrr.nzvals

  @assert nrows == length(rowptrs)-1
  colptrs = Vector{Ti}(undef,ncols+1)
  work = Vector{Ti}(undef,ncols)
  cscnnz = _csrr_to_csc_count!(colptrs,rowptrs,colvals,nzvalscsr,work)
  rowvals = Vector{Ti}(undef,cscnnz)
  nzvalscsc = Vector{Tv}(undef,cscnnz)
  _csrr_to_csc_fill!(colptrs,rowvals,nzvalscsc,rowptrs,colvals,nzvalscsr)
  SparseMatrixCSC(nrows,ncols,colptrs,rowvals,nzvalscsc)
end

# Notation
# csrr: csr with repeated and unsorted columns
# csru: csr witu unsorted columns
# csc: csc with sorted columns

# Adapted form SparseArrays
function _csrr_to_csc_count!(
  colptrs::Vector{Ti},
  rowptrs::Vector{Tj},
  colvals::Vector{Tj},
  nzvalscsr::Vector{Tv},
  work::Vector{Tj}) where {Ti,Tj,Tv}

  nrows = length(rowptrs)-1
  ncols = length(colptrs)-1
  if nrows == 0 || ncols == 0
    fill!(colptrs, Ti(1))
    return Tj(0)
  end

  # Convert csrr to csru by identifying repeated cols with array work.
  # At the same time, count number of unique rows in colptrs shifted by one.
  fill!(colptrs, Ti(0))
  fill!(work, Tj(0))
  writek = Tj(1)
  newcsrrowptri = Ti(1)
  origcsrrowptri = Tj(1)
  origcsrrowptrip1 = rowptrs[2]
  @inbounds for i in 1:nrows
    for readk in origcsrrowptri:(origcsrrowptrip1-Tj(1))
      j = colvals[readk]
      if work[j] < newcsrrowptri
        work[j] = writek
        if writek != readk
          colvals[writek] = j
          nzvalscsr[writek] = nzvalscsr[readk]
        end
        writek += Tj(1)
        colptrs[j+1] += Ti(1)
      else
        klt = work[j]
        nzvalscsr[klt] = +(nzvalscsr[klt], nzvalscsr[readk])
      end
    end
    newcsrrowptri = writek
    origcsrrowptri = origcsrrowptrip1
    origcsrrowptrip1 != writek && (rowptrs[i+1] = writek)
    i < nrows && (origcsrrowptrip1 = rowptrs[i+2])
  end

  # Convert colptrs from counts to ptrs shifted by one
  # (ptrs will be corrected below)
  countsum = Tj(1)
  colptrs[1] = Ti(1)
  @inbounds for j in 2:(ncols+1)
    overwritten = colptrs[j]
    colptrs[j] = countsum
    countsum += overwritten
    @check Base.hastypemax(Ti) && (countsum <= typemax(Ti))
  end

  cscnnz = countsum - Tj(1)
  cscnnz
end

function _csrr_to_csc_fill!(
  colptrs::Vector{Ti},rowvals::Vector{Ti},nzvalscsc::Vector{Tv},
  rowptrs::Vector{Tj},colvals::Vector{Tj},nzvalscsr::Vector{Tv}) where {Ti,Tj,Tv}

  nrows = length(rowptrs)-1
  ncols = length(colptrs)-1
  if nrows == 0 || ncols == 0
    return nothing
  end

  # From csru to csc
  # Tracking write positions in colptrs corrects
  # the column pointers to the final value.
  @inbounds for i in 1:nrows
    for csrk in rowptrs[i]:(rowptrs[i+1]-Tj(1))
      j = colvals[csrk]
      x = nzvalscsr[csrk]
      csck = colptrs[j+1]
      colptrs[j+1] = csck + Ti(1)
      rowvals[csck] = i
      nzvalscsc[csck] = x
    end
  end

  nothing
end


