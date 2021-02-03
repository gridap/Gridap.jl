
function compress(cell_to_bgcell::AbstractVector{<:Integer})
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

function compress(
  cell_to_mat,
  cell_to_bgcell::AbstractVector{<:Integer},
  ccell_to_first_cell::AbstractVector{<:Integer}=compress(cell_to_bgcell))

  nccells = length(ccell_to_first_cell)-1
  ccell_to_cell = lazy_map(Reindex(ccell_to_first_cell),1:nccells)
  ccell_to_bgcell = lazy_map(Reindex(cell_to_bgcell),ccell_to_cell)
  ccell_to_mat = CompressedCellArray(cell_to_mat,ccell_to_first_cell)
  ccell_to_mat, ccell_to_bgcell
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
    for k in eachindex(r)
      r[k] = r[k] + carray[k]
    end
  end
  r
end

@inline function _compress!(matvec::Tuple,cell_to_matvec,cell_range)
  mat, vec = matvec
  rm = copy(mat)
  rv = copy(vec)
  for cell in cell_range
    cmat, cvec = cell_to_matvec[cell]
    for k in eachindex(mat)
      rm[k] = rm[k] + cmat[k]
    end
    for i in eachindex(vec)
      rv[i] = rv[i] + cvec[i]
    end
  end
  (rm,rv)
end

function array_cache(a::CompressedCellArray)
  array = testitem(a.cell_to_array)
  _compress_cache(array,a)
end

function _compress_cache(array,a)
  CachedArray(copy(array)), array_cache(a.cell_to_array)
end

function _compress_cache(matvec::Tuple,a)
  mat, vec = matvec
  cmatvec = (CachedArray(copy(mat)), CachedArray(copy(vec)))
  cmatvec, array_cache(a.cell_to_array)
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
  setsize!(c,size(array))
  r = c.array
  copyto!(r,array)
  for cell in cell_range
    carray = getindex!(ccache,cell_to_array,cell)
    for k in eachindex(array)
      r[k] = r[k] + carray[k]
    end
  end
  r
end

@inline function _compress!(c,matvec::Tuple,ccache,cell_to_matvec,cell_range)
  cm, cv = c
  mat, vec = matvec
  setsize!(cm,size(mat))
  setsize!(cv,size(vec))
  rm = cm.array
  rv = cv.array
  copyto!(rm,mat)
  copyto!(rv,vec)
  for cell in cell_range
    cmat, cvec = getindex!(ccache,cell_to_matvec,cell)
    for k in eachindex(mat)
      rm[k] = rm[k] + cmat[k]
    end
    for i in eachindex(vec)
      rv[i] = rv[i] + cvec[i]
    end
  end
  (rm, rv)
end
