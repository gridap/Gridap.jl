"""
    move_contributions(scell_to_val, strian)
    move_contributions(scell_to_val, strian, ttrian)
    move_contributions(scell_to_val, sglue, tglue)
    move_contributions(scell_to_val, sglue, tglue, nothing)
    move_contributions(scell_to_val, sglue, tglue, mcell_to_scell::AbstractArray)

where `scell_to_val` isa AbstractArray, `s`/`ttrian` isa Triangulation,
`s`/`tglue` isa FaceToFaceGlue.
"""
function move_contributions(
  scell_to_val::AbstractArray, strian::Triangulation)

  scell_to_val, strian
end

function move_contributions(
  scell_to_val::AbstractArray, strian::Triangulation, ttrian::Triangulation)

  D = num_cell_dims(ttrian)
  sglue = get_glue(strian,Val(D))
  tglue = get_glue(ttrian,Val(D))
  move_contributions(scell_to_val,sglue,tglue)
end

function move_contributions(
  scell_to_val::AbstractArray, sglue::FaceToFaceGlue, tglue::FaceToFaceGlue)
  move_contributions(scell_to_val,sglue,tglue,sglue.mface_to_tface)
end

function move_contributions(
  scell_to_val::AbstractArray,
  sglue::FaceToFaceGlue,
  tglue::FaceToFaceGlue,
  mcell_to_scell::AbstractArray)

  mcell_to_val = extend(scell_to_val,mcell_to_scell)
  tcell_to_mcell = tglue.tface_to_mface
  tcell_to_val = lazy_map(Reindex(mcell_to_val),tcell_to_mcell)
end

function move_contributions(
  scell_to_val::AbstractArray,
  sglue::FaceToFaceGlue,
  tglue::FaceToFaceGlue,
  mcell_to_scell::Nothing)

  scell_to_mcell = sglue.tface_to_mface
  mcell_to_tcell = tglue.mface_to_tface
  ntcells = length(tglue.tface_to_mface)
  tcell_to_scells_ptr = zeros(Int32,ntcells+1)
  for mcell in scell_to_mcell
    tcell = mcell_to_tcell[mcell]
    tcell_to_scells_ptr[1+tcell] += Int32(1)
  end
  length_to_ptrs!(tcell_to_scells_ptr)
  ndata = tcell_to_scells_ptr[end]-1
  tcell_to_scells_data = zeros(Int32,ndata)
  for (scell,mcell) in enumerate(scell_to_mcell)
    tcell = mcell_to_tcell[mcell]
    p = tcell_to_scells_ptr[tcell]
    tcell_to_scells_data[p] = scell
    tcell_to_scells_ptr[tcell] += Int32(1)
  end
  rewind_ptrs!(tcell_to_scells_ptr)
  tcell_to_scells = Table(tcell_to_scells_data,tcell_to_scells_ptr)
  k = CombineContributionsMap(scell_to_val)
  tcell_to_val = lazy_map(k,tcell_to_scells)
end

struct CombineContributionsMap{A} <: Map
  scell_to_val::A
end

function return_cache(k::CombineContributionsMap{<:AbstractArray{<:Number}},scells)
  z = zero(testitem(k.scell_to_val))
  z, array_cache(k.scell_to_val)
end

function evaluate!(cache,k::CombineContributionsMap{<:AbstractArray{<:Number}},scells)
  z,c = cache
  val = zero(z)
  for scell in scells
    val += getindex!(c,k.scell_to_val,scell)
  end
  val
end

function return_cache(k::CombineContributionsMap,scells)
  array = testitem(k.scell_to_val)
  c12 = _cache_compress(array)
  c = array_cache(k.scell_to_val)
  c12, c
end

function evaluate!(cache,k::CombineContributionsMap,scells)
  c12,c = cache
  c1,c2 = c12
  if length(scells) == 0
    _setempty_compress!(c1)
    val = _uncached_compress!(c1,c2)
  else
    val1 = getindex!(c,k.scell_to_val,scells[1])
    _setsize_compress!(c1,val1)
    val = _uncached_compress!(c1,c2)
    _copyto_compress!(val,val1)
    for i in 2:length(scells)
      vali = getindex!(c,k.scell_to_val,scells[i])
      _addto_compress!(val,vali)
    end
    val
  end
  val
end

function _cache_compress(array::AbstractArray)
  c1 = CachedArray(copy(array))
  c1, nothing
end

function _cache_compress(array::ArrayBlock)
  c1 = CachedArray(deepcopy(array))
  c2 = return_cache(Fields.unwrap_cached_array,c1)
  c1,c2
end

function _cache_compress(matvec::Tuple)
  mat, vec = matvec
  cmat,ca = _cache_compress(mat)
  cvec,cb = _cache_compress(vec)
  (cmat, cvec), (ca,cb)
end

function _setempty_compress!(a::CachedArray)
  setsize!(a,0 .* size(a))
end

function _setempty_compress!(a::ArrayBlock)
  for k in eachindex(a.array)
    if a.touched[k]
      _setempty_compress!(a.array[k])
    end
  end
end

function _setempty_compress!(ab::Tuple)
  a,b = ab
  _setempty_compress!(a)
  _setempty_compress!(b)
end

function _uncached_compress!(a::CachedArray,c2)
  a.array
end

function _uncached_compress!(c1::ArrayBlock,c2)
  evaluate!(c2,Fields.unwrap_cached_array,c1)
end

function _uncached_compress!(c1::Tuple,c2)
  a1,b1 = c1
  a2,b2 = c2
  a = _uncached_compress!(a1,a2)
  b = _uncached_compress!(b1,b2)
  (a,b)
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

function _setsize_compress!(a::Tuple,b::Tuple)
  a1,a2 = a
  b1,b2 = b
  _setsize_compress!(a1,b1)
  _setsize_compress!(a2,b2)
end

function _copyto_compress!(a::AbstractArray,b::AbstractArray)
  @check size(a) == size(b)
  for k in eachindex(a)
    a[k] = b[k]
  end
  a
end

function _copyto_compress!(a::ArrayBlock,b::ArrayBlock)
  @check size(a) == size(b)
  for k in eachindex(a.array)
    @check a.touched[k] == b.touched[k]
    if a.touched[k]
      _copyto_compress!(a.array[k],b.array[k])
    end
  end
  a
end

function _copyto_compress!(a::Tuple,b::Tuple)
  a1,a2 = a
  b1,b2 = b
  _copyto_compress!(a1,b1)
  _copyto_compress!(a2,b2)
end

function _addto_compress!(a::AbstractArray,b::AbstractArray)
  @check size(a) == size(b)
  for k in eachindex(a)
    a[k] = a[k] + b[k]
  end
  a
end

function _addto_compress!(a::ArrayBlock,b::ArrayBlock)
  @check size(a) == size(b)
  for k in eachindex(a.array)
    @check a.touched[k] == b.touched[k]
    if a.touched[k]
      _addto_compress!(a.array[k],b.array[k])
    end
  end
  a
end

function _addto_compress!(a::Tuple,b::Tuple)
  a1,a2 = a
  b1,b2 = b
  _addto_compress!(a1,b1)
  _addto_compress!(a2,b2)
end

