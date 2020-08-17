
function attach_dirichlet(cellmatvec,cellvals,cellmask=Fill(true,length(cellvals)))
  k = AttachDirichletKernel()
  apply(k,cellmatvec,cellvals,cellmask)
end

struct AttachDirichletKernel <: Kernel
  muladd::MulAddKernel{Int}
  AttachDirichletKernel() = new(MulAddKernel(-1,1))
end

function Arrays.kernel_cache(k::AttachDirichletKernel,matvec::Tuple,vals,mask)
  mat, vec = matvec
  kernel_cache(k.muladd,mat,vals,vec)
end

@inline function Arrays.apply_kernel!(cache,k::AttachDirichletKernel,matvec::Tuple,vals,mask)
  if mask
    mat, vec = matvec
    vec_with_bcs = apply_kernel!(cache,k.muladd,mat,vals,vec)
    (mat, vec_with_bcs)
  else
    matvec
  end
end
