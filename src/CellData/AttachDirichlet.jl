
function attach_dirichlet(cellmatvec,cellvals,cellmask=Fill(true,length(cellvals)))
  k = AttachDirichletMap()
  lazy_map(k,cellmatvec,cellvals,cellmask)
end

struct AttachDirichletMap <: Map
  muladd::MulAddMap{Int}
  AttachDirichletMap() = new(MulAddMap(-1,1))
end

function Arrays.return_cache(k::AttachDirichletMap,matvec::Tuple,vals,mask)
  mat, vec = matvec
  return_cache(k.muladd,mat,vals,vec)
end

@inline function Arrays.evaluate!(cache,k::AttachDirichletMap,matvec::Tuple,vals,mask)
  if mask
    mat, vec = matvec
    vec_with_bcs = evaluate!(cache,k.muladd,mat,vals,vec)
    (mat, vec_with_bcs)
  else
    matvec
  end
end

function Arrays.return_cache(k::AttachDirichletMap,mat::AbstractMatrix,vals,mask)
  cm = return_cache(*,mat,vals)
  cv = CachedArray(mat*vals)
  fill!(cv.array,zero(eltype(cv)))
  (cm,cv)
end

@inline function Arrays.evaluate!(cache,k::AttachDirichletMap,mat::AbstractMatrix,vals,mask)
  cm, cv = cache
  if mask
    vec_with_bcs = evaluate!(cm,*,mat,vals)
    scale_entries!(vec_with_bcs,-1)
    (mat, vec_with_bcs)
  else
    if size(mat,1) != size(cv,1)
      m = axes(mat,1)
      setaxes!(cv,(m,))
      fill!(cv.array,zero(eltype(cv)))
    end
    (mat, cv.array)
  end
end
