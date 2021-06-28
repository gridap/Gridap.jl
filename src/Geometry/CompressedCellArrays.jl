
function compress_ids(cell_to_bgcell::AbstractVector{<:Integer})
  nccells = 0
  c = -1
  for bgcell in cell_to_bgcell
    if c != bgcell
      nccells += 1
      c = bgcell
    end
  end
  ccell_to_first_cell  = zeros(Int32,nccells+1)
  ccell = 0
  c = -1
  for (cell,bgcell) in enumerate(cell_to_bgcell)
    if c != bgcell
      ccell += 1
      c = bgcell
      ccell_to_first_cell[ccell] = cell
    end
  end
  ccell_to_first_cell[end] = length(cell_to_bgcell)+1
  ccell_to_first_cell
end

function compress_ids(cell_ids,cell_to_bgcell::AbstractVector{<:Integer})
  ccell_to_first_cell = compress_ids(cell_to_bgcell)
  nccells = length(ccell_to_first_cell)-1
  ccell_to_cell = view(ccell_to_first_cell,1:nccells)
  ccell_ids = view(cell_ids,ccell_to_cell)
  ccell_ids
end

function compress_contributions(
  cell_to_mat,
  cell_to_bgcell::AbstractVector{<:Integer},
  ccell_to_first_cell::AbstractVector{<:Integer}=compress_ids(cell_to_bgcell))

  nccells = length(ccell_to_first_cell)-1
  ccell_to_cell = view(ccell_to_first_cell,1:nccells)
  ccell_to_mat = CompressedCellArray(cell_to_mat,ccell_to_first_cell)
  ccell_to_mat
end

struct CompressedCellArray{T,A,B} <: AbstractVector{T}
  cell_to_array::A
  ccell_to_first_cell::B
  function CompressedCellArray(cell_to_array::A, ccell_to_first_cell::B) where {A,B}
    T = eltype(cell_to_array)
    new{T,A,B}(cell_to_array,ccell_to_first_cell)
  end
end

Base.size(a::CompressedCellArray) = (length(a.ccell_to_first_cell)-1,)
Base.IndexStyle(::Type{<:CompressedCellArray}) = IndexLinear()

function Base.getindex(a::CompressedCellArray,ccell::Integer)
  cell_first = a.ccell_to_first_cell[ccell]
  cell_last = a.ccell_to_first_cell[ccell+1]-1
  array = a.cell_to_array[cell_first]
  cell_range = (cell_first+1):cell_last
  _compress!(array,a.cell_to_array,cell_range)
end

@inline function _compress!(array,cell_to_array,cell_range)
  r = copy(array)
  for cell in cell_range
    carray = cell_to_array[cell]
    _addto!(r,carray)
  end
  r
end

@inline function _compress!(matvec::Tuple,cell_to_matvec,cell_range)
  mat, vec = matvec
  rm = copy(mat)
  rv = copy(vec)
  for cell in cell_range
    cmat, cvec = cell_to_matvec[cell]
    _addto!(rm,cmat)
    _addto!(rv,cvec)
  end
  (rm,rv)
end

function _addto!(a,b)
  @check size(a) == size(b)
  for k in eachindex(a)
    a[k] = a[k] + b[k]
  end
  a
end

function _addto!(a::ArrayBlock,b::ArrayBlock)
  @check size(a) == size(b)
  for k in eachindex(a.array)
    @check a.touched[k] == b.touched[k]
    if a.touched[k]
      _addto!(a.array[k],b.array[k])
    end
  end
  a
end

function array_cache(a::CompressedCellArray)
  array = testitem(a.cell_to_array)
  c = _compress_cache(array)
  ccache = array_cache(a.cell_to_array)
  c,ccache
end

function _compress_cache(array)
  c1 = CachedArray(copy(array))
  c2 = return_cache(Fields.unwrap_cached_array,c1)
  c1,c2
end

function _compress_cache(matvec::Tuple)
  mat, vec = matvec
  cmat = _compress_cache(mat)
  cvec = _compress_cache(vec)
  cmat, cvec
end

function getindex!(cache,a::CompressedCellArray,ccell::Integer)
  c, ccache = cache
  cell_first = a.ccell_to_first_cell[ccell]
  cell_last = a.ccell_to_first_cell[ccell+1]-1
  array = getindex!(ccache,a.cell_to_array,cell_first)
  cell_range = (cell_first+1):cell_last
  _compress!(c,array,ccache,a.cell_to_array,cell_range)
end

@inline function _compress!(c,array,ccache,cell_to_array,cell_range)
  c1,c2 = c
  _setsize_compress!(c1,array)
  r = evaluate!(c2,Fields.unwrap_cached_array,c1)
  copyto!(r,array)
  for cell in cell_range
    carray = getindex!(ccache,cell_to_array,cell)
    _addto!(r,carray)
  end
  r
end

@inline function _compress!(c,matvec::Tuple,ccache,cell_to_matvec,cell_range)
  cm, cv = c
  cm1,cm2 = cm
  cv1,cv2 = cv
  mat, vec = matvec
  _setsize_compress!(cm1,mat)
  _setsize_compress!(cv1,vec)
  rm = evaluate!(cm2,Fields.unwrap_cached_array,cm1)
  rv = evaluate!(cv2,Fields.unwrap_cached_array,cv1)
  copyto!(rm,mat)
  copyto!(rv,vec)
  for cell in cell_range
    cmat, cvec = getindex!(ccache,cell_to_matvec,cell)
    _addto!(rm,cmat)
    _addto!(rv,cvec)
  end
  (rm, rv)
end

function _setsize_compress!(a::CachedArray,b::AbstractArray)
  setsize!(a,size(b))
end

function _setsize_compress!(a::ArrayBlock,b::ArrayBlock)
  @check size(a) == size(b)
  for k in eachindex(a.array)
    @check a.touched[k] == b.touched[k]
    if a.touched[k]
    _setsize_compress!(a.array[k],b.array[k])
    end
  end
end


