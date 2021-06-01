
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

struct CounterCSC{Tv,Ti,L}
  tv::Type{Tv}
  nrows::Int
  ncols::Int
  colnnzmax::Vector{Ti}
  loop_style::L
end

LoopStyle(::Type{CounterCSC{Tv,Ti,L}}) where {Tv,Ti,L} = L()

@inline function add_entry!(::typeof(+),a::CounterCSC{Tv,Ti,Loop},v,i,j) where {Tv,Ti}
  a.colnnzmax[j] += Ti(1)
  nothing
end

@inline function add_entry!(::typeof(+),a::CounterCSC{Tv,Ti,DoNotLoop},v,i,j) where {Tv,Ti}
  nothing
end

struct InserterCSC{Tv,Ti}
  nrows::Int
  ncols::Int
  colptr::Vector{Ti}
  colnnz::Vector{Ti}
  rowval::Vector{Ti}
  nzval::Vector{Tv}
end

LoopStyle(::Type{<:InserterCSC}) = Loop()

@inline function add_entry!(::typeof(+),a::InserterCSC{Tv,Ti},v::Nothing,i,j) where {Tv,Ti}
  pini = Int(a.colptr[j])
  pend = pini + Int(a.colnnz[j]) - 1
  p = searchsortedfirst(a.rowval,i,pini,pend,Base.Order.Forward)
  if (p>pend)
    # add new entry
    a.colnnz[j] += 1
    a.rowval[p] = i
  elseif a.rowval[p] != i
    # shift one forward from p to pend
    @check  pend+1 < Int(a.colptr[j+1])
    for k in pend:-1:p
      o = k + 1
      a.rowval[o] = a.rowval[k]
    end
    # add new entry
    a.colnnz[j] += 1
    a.rowval[p] = i
  end
  nothing
end

@noinline function add_entry!(::typeof(+),a::InserterCSC{Tv,Ti},v,i,j) where {Tv,Ti}
  pini = Int(a.colptr[j])
  pend = pini + Int(a.colnnz[j]) - 1
  p = searchsortedfirst(a.rowval,i,pini,pend,Base.Order.Forward)
  if (p>pend)
    # add new entry
    a.colnnz[j] += 1
    a.rowval[p] = i
    a.nzval[p] = v
  elseif a.rowval[p] != i
    # shift one forward from p to pend
    @check  pend+1 < Int(a.colptr[j+1])
    for k in pend:-1:p
      o = k + 1
      a.rowval[o] = a.rowval[k]
      a.nzval[o] = a.nzval[k]
    end
    # add new entry
    a.colnnz[j] += 1
    a.rowval[p] = i
    a.nzval[p] = v
  else 
    # update existing entry
    a.nzval[p] += v 
  end
  nothing
end

## index of the first value of vector a that is greater than or equal to x;
## returns length(v)+1 if x is greater than all values in v.
#function _searchsortedfirst(v, x, lo::T, hi::T) where T<:Integer
#  if x <= v[lo]
#    return lo
#  end
#  u = one(T)
#  if x > v[hi]
#    return hi + u
#  end
#  d = T(10) - u
#  @inbounds while lo < hi - d
#    m = Base.Sort.midpoint(lo, hi)
#    if v[m] < x
#      lo = m
#    else
#      hi = m
#    end
#  end
#  i = lo
#  if x<=v[i]
#    return i
#  end
#  i += u
#  if x<=v[i]
#    return i
#  end
#  i += u
#  if x<=v[i]
#    return i
#  end
#  i += u
#  if x<=v[i]
#    return i
#  end
#  i += u
#  if x<=v[i]
#    return i
#  end
#  i += u
#  if x<=v[i]
#    return i
#  end
#  i += u
#  if x<=v[i]
#    return i
#  end
#  i += u
#  if x<=v[i]
#    return i
#  end
#  i += u
#  if x<=v[i]
#    return i
#  end
#  i += u
#  if x<=v[i]
#    return i
#  end
#  return hi + u
#end

function nz_counter(
  builder::SparseMatrixBuilder{SparseMatrixCSC{Tv,Ti},<:MinMemory},
  axes) where {Tv,Ti}
  nrows = length(axes[1])
  ncols = length(axes[2])
  maxnnz = builder.approach.maxnnz
  if isa(maxnnz,Nothing)
    colnnzmax = zeros(Ti,ncols)
    CounterCSC(Tv,nrows,ncols,colnnzmax,Loop())
  elseif isa(maxnnz,Integer)
    colnnzmax = fill(Ti(maxnnz),ncols)
    CounterCSC(Tv,nrows,ncols,colnnzmax,DoNotLoop())
  elseif isa(maxnnz,Vector{<:Integer})
    colnnzmax = maxnnz
    @assert length(colnnzmax) == ncols
    CounterCSC(Tv,nrows,ncols,colnnzmax,DoNotLoop())
  else
    @notimplemented
  end
end

function nz_counter(
  builder::SparseMatrixBuilder{SparseMatrixCSC{Tv,Ti},<:MinCPU},
  axes) where {Tv,Ti}

  nrows = length(axes[1])
  ncols = length(axes[2])
  rowptrs = zeros(Ti,nrows+1)
  CounterCSRR(Tv,nrows,ncols,rowptrs)
end

#function nz_counter(::Type{SparseMatrixCSC{Tv,Ti}},axes) where {Tv,Ti}
#  builder = SparseMatrixBuilder(SparseMatrixCSC{Tv,Ti},MinMemory())
#  nz_counter(builder,axes)
#end

function nz_allocation(a::CounterCSC{Tv,Ti}) where {Tv,Ti}
  colptr = Vector{Ti}(undef,a.ncols+1)
  @inbounds for i in 1:a.ncols
    colptr[i+1] = a.colnnzmax[i]
  end
  length_to_ptrs!(colptr)
  ndata = colptr[end] - one(Ti)
  rowval = Vector{Ti}(undef,ndata)
  nzval = zeros(Tv,ndata)
  colnnz = a.colnnzmax
  fill!(colnnz,zero(Ti))
  InserterCSC(a.nrows,a.ncols,colptr,colnnz,rowval,nzval)
end

function create_from_nz(a::InserterCSC{Tv,Ti}) where {Tv,Ti}
  k = 1
  for j in 1:a.ncols
    pini = Int(a.colptr[j])
    pend = pini + Int(a.colnnz[j]) - 1
    for p in pini:pend
      a.nzval[k] = a.nzval[p]
      a.rowval[k] = a.rowval[p]
      k += 1
    end
  end
  @inbounds for j in 1:a.ncols
    a.colptr[j+1] = a.colnnz[j]
  end
  length_to_ptrs!(a.colptr)
  nnz = a.colptr[end]-1
  resize!(a.rowval,nnz)
  resize!(a.nzval,nnz)
  SparseMatrixCSC(a.nrows,a.ncols,a.colptr,a.rowval,a.nzval)
end

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
  work::Vector{Ti}
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

function nz_allocation(a::CounterCSRR{Tv,Ti}) where {Tv,Ti}
  rowptrs = a.rowptrs
  length_to_ptrs!(rowptrs)
  ndata = rowptrs[end]-1
  colvals = Vector{Ti}(undef,ndata)
  nzvals = zeros(Tv,ndata)
  work = Vector{Ti}(undef,a.ncols)
  CSRR(a.nrows,a.ncols,rowptrs,colvals,nzvals,work)
end

function create_from_nz(a::CSRR{Tv,Ti}) where {Tv,Ti}
  rewind_ptrs!(a.rowptrs)
  colptrs = Vector{Ti}(undef,a.ncols+1)
  cscnnz = _csrr_to_csc_count!(colptrs,a.rowptrs,a.colvals,a.nzvals,a.work)
  rowvals = Vector{Ti}(undef,cscnnz)
  nzvalscsc = Vector{Tv}(undef,cscnnz)
  _csrr_to_csc_fill!(colptrs,rowvals,nzvalscsc,a.rowptrs,a.colvals,a.nzvals)
  SparseMatrixCSC(a.nrows,a.ncols,colptrs,rowvals,nzvalscsc)
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


