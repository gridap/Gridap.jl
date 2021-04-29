
function attach_dirichlet(cellmatvec,cellvals,cellmask=Fill(true,length(cellvals)))
  k = AttachDirichletMap()
  lazy_map(k,cellmatvec,cellvals,cellmask)
end

struct ZeroVectorMap <: Map end

function Arrays.return_cache(k::ZeroVectorMap,mat::AbstractMatrix)
  v = similar(mat,eltype(mat),size(mat,1))
  fill!(v,zero(eltype(v)))
  CachedArray(v)
end

function Arrays.evaluate!(cache,k::ZeroVectorMap,mat::AbstractMatrix)
  setsize!(cache,(size(mat,1),))
  v = cache.array
  fill!(v,zero(eltype(v)))
  v
end

function Arrays.return_cache(k::ZeroVectorMap,mat::MatrixBlock)
  ni,nj = size(mat)
  ai = testitem(mat)
  ci = return_cache(k,ai)
  vi = evaluate!(ci,k,ai)
  c = Vector{typeof(ci)}(undef,ni)
  array = Vector{typeof(vi)}(undef,ni)
  touched = fill(false,ni)
  for i in 1:ni
    for j in 1:nj
      if !touched[i] && mat.touched[i,j]
        c[i] = return_cache(k,mat.array[i,j])
        touched[i] = true
      end
    end
  end
  ArrayBlock(array,touched), c
end

function Arrays.evaluate!(cache,k::ZeroVectorMap,mat::MatrixBlock)
  r,c = cache
  ni,nj = size(mat)
  fill!(r.touched,false)
  for j in 1:nj
    for i in 1:ni
      if !r.touched[i] && mat.touched[i,j]
        r.array[i] = evaluate!(c[i],k,mat.array[i,j])
        r.touched[i] = true
      end
    end
  end
  r
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

function Arrays.return_cache(k::AttachDirichletMap,mat,vals,mask)
  cm = return_cache(*,mat,vals)
  cv = return_cache(ZeroVectorMap(),mat)
  (cm,cv)
end

@inline function Arrays.evaluate!(cache,k::AttachDirichletMap,mat,vals,mask)
  cm, cv = cache
  if mask
    vec_with_bcs = evaluate!(cm,*,mat,vals)
    scale_entries!(vec_with_bcs,-1)
  else
    vec_with_bcs = evaluate!(cv,ZeroVectorMap(),mat)
  end
  (mat, vec_with_bcs)
end


