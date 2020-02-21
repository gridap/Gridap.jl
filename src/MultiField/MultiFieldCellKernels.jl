
function kernel_cache(k::CellMatVecKernel,a::MultiFieldArray,b::MultiFieldArray,c...)
  mat = _prepare_multifield_mat(a,b)
  vec = _prepare_multifield_vec(a)
  cmat = CachedMultiFieldArray(mat)
  cvec = CachedMultiFieldArray(vec)
  (mat, vec, cmat, cvec)
end

@inline function apply_kernel!(cache,k::CellMatVecKernel,a::MultiFieldArray,b::MultiFieldArray,c...)
  mat, vec, cmat, cvec = cache
  _resize_multifield_mat!(cmat,a,b)
  _resize_multifield_vec!(cvec,a)
  _move_cached_arrays!(mat,cmat)
  _move_cached_arrays!(vec,cvec)
  z = zero(eltype(mat))
  fill!(mat,z)
  fill!(vec,z)
  k.op(mat,vec,a,b,c...)
  (mat,vec)
end

function kernel_cache(k::CellMatKernel,a::MultiFieldArray,b::MultiFieldArray,c...)
  mat = _prepare_multifield_mat(a,b)
  cmat = CachedMultiFieldArray(mat)
  (mat, cmat)
end

@inline function apply_kernel!(cache,k::CellMatKernel,a::MultiFieldArray,b::MultiFieldArray,c...)
  mat, cmat = cache
  _resize_multifield_mat!(cmat,a,b)
  _move_cached_arrays!(mat,cmat)
  z = zero(eltype(mat))
  fill!(mat,z)
  k.op(mat,a,b,c...)
  mat
end

function kernel_cache(k::CellVecKernel,a::MultiFieldArray,c...)
  vec = _prepare_multifield_vec(a)
  cvec = CachedMultiFieldArray(vec)
  (vec, cvec)
end

@inline function apply_kernel!(cache,k::CellVecKernel,a::MultiFieldArray,c...)
  vec, cvec = cache
  _resize_multifield_vec!(cvec,a)
  _move_cached_arrays!(vec,cvec)
  z = zero(eltype(vec))
  fill!(vec,z)
  k.op(vec,a,c...)
  vec
end

function _prepare_multifield_mat(a,b)

  function fun1(i,j,a,b)
    ai = a.blocks[i]
    bj = b.blocks[j]
    Ta = eltype(ai)
    Tb = eltype(bj)
    T = return_type(inner,Ta,Tb)
    m = size(ai,2)
    n = size(bj,2)
    zeros(T,m,n)
  end

  function fun2(i,j,a,b)
    ci, = a.coordinates[i]
    cj, = b.coordinates[j]
    (ci,cj)
  end

  na = length(a.blocks)
  nb = length(b.blocks)
  blocks = [ fun1(i,j,a,b) for i in 1:na for j in 1:nb ]
  coordinates = [ fun2(i,j,a,b) for i in 1:na for j in 1:nb ]
  MultiFieldArray(blocks,coordinates)
end

function _prepare_multifield_vec(a)

  function fun1(i,a)
    ai = a.blocks[i]
    Ta = eltype(ai)
    T = return_type(inner,Ta,Ta)
    m = size(ai,2)
    zeros(T,m)
  end

  function fun2(i,a)
    ci, = a.coordinates[i]
    (ci,)
  end

  na = length(a.blocks)
  blocks = [ fun1(i,a) for i in 1:na ]
  coordinates = [ fun2(i,a) for i in 1:na ]
  MultiFieldArray(blocks,coordinates)
end

function _resize_multifield_mat!(c,a,b)
  for i in 1:length(a.blocks)
    for j in 1:length(b.blocks)
      ai = a.blocks[i]
      bj = b.blocks[j]
      m = size(ai,2)
      n = size(bj,2)
      ci, = a.coordinates[i]
      cj, = b.coordinates[j]
      k = c.ptrs[ci,cj]
      ck = c.blocks[k]
      setsize!(ck,(m,n))
    end
  end
end

function _resize_multifield_vec!(c,a)
  for i in 1:length(a.blocks)
    ai = a.blocks[i]
    m = size(ai,2)
    ci, = a.coordinates[i]
    k = c.ptrs[ci]
    ck = c.blocks[k]
    setsize!(ck,(m,))
  end
end
