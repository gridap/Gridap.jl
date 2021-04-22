
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
    _zero_if_needed!(cv,mat)
    (mat, cv.array)
  end
end

function Arrays.return_cache(k::AttachDirichletMap,mat::GBlock,vals,mask)
  cm = return_cache(*,mat,vals)
  ai = testvalue(eltype(mat))
  vi = testvalue(eltype(vals))
  _, cvi = return_cache(k,ai,vi,mask)
  array = Vector{typeof(cvi)}(undef,length(cm))
  touched = fill(false,length(cm))
  for i in 1:size(mat.array,1)
    for j in 1:size(mat.array,2)
      if mat.touched[i,j] && vals.touched[j]
        if !touched[i]
          _, cvi = return_cache(k,mat.array[i,j],vals.array[j],mask)
          array[i] = cvi
          touched[i] = true
        end
      end
    end
  end
  (cm,GBlock(array,touched))
end

@inline function Arrays.evaluate!(cache,k::AttachDirichletMap,mat::GBlock,vals,mask)
  cm, cv = cache
  if mask
    vec_with_bcs = evaluate!(cm,*,mat,vals)
    scale_entries!(vec_with_bcs,-1)
    (mat, vec_with_bcs)
  else
    _zero_if_needed!(cv,mat)
    (mat, cv)
  end
end

function _zero_if_needed!(cv::CachedVector,mat::AbstractMatrix)
  if size(mat,1) != size(cv,1)
    m = axes(mat,1)
    setaxes!(cv,(m,))
    fill!(cv.array,zero(eltype(cv)))
  end
  nothing
end

function _zero_if_needed!(cv::GBlock,mat::GBlock)
  ni, nj = size(mat.array)
  for i in 1:ni
    for j in 1:nj
      if mat.touched[i,j]
        _zero_if_needed!(cv.array[i],mat.array[i,j])
        break
      end
    end
  end
  nothing
end


