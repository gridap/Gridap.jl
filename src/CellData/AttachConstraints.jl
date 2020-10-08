
function attach_constraints_rows(cellvec,cellconstr,cellmask=Fill(true,length(cellconstr)))
  lazy_map(ConstrainRowsMap(),cellvec,cellconstr,cellmask)
end

function attach_constraints_cols(cellmat,cellconstr,cellmask=Fill(true,length(cellconstr)))
  cellconstr_t = lazy_map(transpose,cellconstr)
  lazy_map(ConstrainColsMap(),cellmat,cellconstr_t,cellmask)
end

struct ConstrainRowsMap <: Map end

function Arrays.return_cache(k::ConstrainRowsMap,array::AbstractArray,constr,mask)
  return_cache(MulMap(),constr,array)
end

@inline function Arrays.evaluate!(cache,k::ConstrainRowsMap,array::AbstractArray,constr,mask)
  if mask
    evaluate!(cache,MulMap(),constr,array)
  else
    array
  end
end

function Arrays.return_cache(k::ConstrainRowsMap,matvec::Tuple,constr,mask)
  mat, vec = matvec
  cmat = return_cache(k,mat,constr,mask)
  cvec = return_cache(k,vec,constr,mask)
  (cmat,cvec)
end

@inline function Arrays.evaluate!(cache,k::ConstrainRowsMap,matvec::Tuple,constr,mask)
  if mask
    cmat, cvec = cache
    mat, vec = matvec
    _mat = evaluate!(cmat,k,mat,constr,mask)
    _vec = evaluate!(cvec,k,vec,constr,mask)
    (_mat,_vec)
  else
    matvec
  end
end

struct ConstrainColsMap <: Map end

function Arrays.return_cache(k::ConstrainColsMap,array::AbstractArray,constr_t,mask)
  return_cache(MulMap(),array,constr_t)
end

@inline function Arrays.evaluate!(cache,k::ConstrainColsMap,array::AbstractArray,constr_t,mask)
  if mask
    evaluate!(cache,MulMap(),array,constr_t)
  else
    array
  end
end

function Arrays.return_cache(k::ConstrainColsMap,matvec::Tuple,constr_t,mask)
  mat, vec = matvec
  return_cache(k,mat,constr_t,mask)
end

@inline function Arrays.evaluate!(cache,k::ConstrainColsMap,matvec::Tuple,constr_t,mask)
  if mask
    mat, vec = matvec
    _mat = evaluate!(cache,k,mat,constr_t,mask)
    (_mat,vec)
  else
    matvec
  end
end

function merge_cell_constraints_at_skeleton(cL,cR,axesL_rows,axesR_rows,axesL_cols,axesR_cols)
  blocks = (cL,cR)
  blockids = [(1,1),(2,2)]
  axs_rows = create_array_of_blocked_axes(axesL_rows,axesR_rows)
  axs_cols = create_array_of_blocked_axes(axesL_cols,axesR_cols)
  axs = lazy_map((r,c) -> (r[1],c[1]),axs_rows,axs_cols)
  VectorOfBlockArrayCoo(blocks,blockids,axs)
end

function identity_constraints(cell_axes)
  lazy_map(IdentityConstraintMap(),cell_axes)
end

struct IdentityConstraintMap <: Map end

function Arrays.return_cache(k::IdentityConstraintMap,axs)
  n = length(axs[1])
  a = zeros(n,n)
  CachedArray(a)
end

function Arrays.evaluate!(cache,k::IdentityConstraintMap,axs)
  n = length(axs[1])
  setsize!(cache,(n,n))
  a = cache.array
  fill!(a,zero(eltype(a)))
  o = one(eltype(a))
  @inbounds for i in 1:size(a,1)
    a[i,i] = o
  end
  a
end
