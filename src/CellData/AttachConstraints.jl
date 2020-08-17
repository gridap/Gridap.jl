
function attach_constraints_rows(cellvec,cellconstr,cellmask=Fill(true,length(cellconstr)))
  apply(ConstrainRowsKernel(),cellvec,cellconstr,cellmask)
end

function attach_constraints_cols(cellmat,cellconstr,cellmask=Fill(true,length(cellconstr)))
  cellconstr_t = apply(transpose,cellconstr)
  apply(ConstrainColsKernel(),cellmat,cellconstr_t,cellmask)
end

struct ConstrainRowsKernel <: Kernel end

function Arrays.kernel_cache(k::ConstrainRowsKernel,array::AbstractArray,constr,mask)
  kernel_cache(MulKernel(),constr,array)
end

@inline function Arrays.apply_kernel!(cache,k::ConstrainRowsKernel,array::AbstractArray,constr,mask)
  if mask
    apply_kernel!(cache,MulKernel(),constr,array)
  else
    array
  end
end

function Arrays.kernel_cache(k::ConstrainRowsKernel,matvec::Tuple,constr,mask)
  mat, vec = matvec
  cmat = kernel_cache(k,mat,constr,mask)
  cvec = kernel_cache(k,vec,constr,mask)
  (cmat,cvec)
end

@inline function Arrays.apply_kernel!(cache,k::ConstrainRowsKernel,matvec::Tuple,constr,mask)
  if mask
    cmat, cvec = cache
    mat, vec = matvec
    _mat = apply_kernel!(cmat,k,mat,constr,mask)
    _vec = apply_kernel!(cvec,k,vec,constr,mask)
    (_mat,_vec)
  else
    matvec
  end
end

struct ConstrainColsKernel <: Kernel end

function Arrays.kernel_cache(k::ConstrainColsKernel,array::AbstractArray,constr_t,mask)
  kernel_cache(MulKernel(),array,constr_t)
end

@inline function Arrays.apply_kernel!(cache,k::ConstrainColsKernel,array::AbstractArray,constr_t,mask)
  if mask
    apply_kernel!(cache,MulKernel(),array,constr_t)
  else
    array
  end
end

function Arrays.kernel_cache(k::ConstrainColsKernel,matvec::Tuple,constr_t,mask)
  mat, vec = matvec
  kernel_cache(k,mat,constr_t,mask)
end

@inline function Arrays.apply_kernel!(cache,k::ConstrainColsKernel,matvec::Tuple,constr_t,mask)
  if mask
    mat, vec = matvec
    _mat = apply_kernel!(cache,k,mat,constr_t,mask)
    (_mat,vec)
  else
    matvec
  end
end

function merge_cell_constraints_at_skeleton(cL,cR,axesL_rows,axesR_rows,axesL_cols,axesR_cols)
  blocks = (cL,cR)
  blockids = [(1,1),(2,2)]
  axs_rows = apply(Fields._cat_axes,axesL_rows,axesR_rows)
  axs_cols = apply(Fields._cat_axes,axesL_cols,axesR_cols)
  axs = apply((r,c) -> (r[1],c[1]),axs_rows,axs_cols)
  VectorOfBlockArrayCoo(blocks,blockids,axs)
end

