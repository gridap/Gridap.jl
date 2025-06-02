
# Make Base.:* behave like a Map

return_value(::typeof(*),a::Number,b::AbstractArray{<:Number}) = a*b

function return_cache(::typeof(*),a::Number,b::AbstractArray{<:Number})
  T = typeof(a*testitem(b))
  N = ndims(b)
  CachedArray(T,N)
end

function evaluate!(cache,::typeof(*),a::Number,b::AbstractArray{<:Number})
  setsize!(cache,size(b))
  c = cache.array
  rmul!(b,a)
  return c
end

return_value(::typeof(*),a::AbstractArray{<:Number},b::Number) = a*b
return_cache(::typeof(*),a::AbstractArray{<:Number},b::Number) = return_cache(*,b,a)
evaluate!(cache,::typeof(*),a::AbstractArray{<:Number},b::Number) = evaluate!(cache,*,b,a)

function return_cache(::typeof(*),a::AbstractArray{<:Number},b::AbstractArray{<:Number})
  T = typeof(testitem(a)*testitem(b))
  N = ifelse(isa(b,AbstractVector),1,2)
  CachedArray(T,N)
end

function evaluate!(cache,::typeof(*),a::AbstractArray{<:Number},b::AbstractArray{<:Number})
  setsize_op!(*,cache,a,b)
  c = cache.array
  mul!(c,a,b)
  return c
end

function return_value(::typeof(*),a::ArrayBlock{A,2},b::ArrayBlock{B,1}) where {A,B}
  ri = return_value(*,testvalue(A),testvalue(B))
  array = Vector{typeof(ri)}(undef,size(a.array,1))
  touched = fill(false,size(a.array,1))
  ArrayBlock(array,touched)
end

function return_cache(::typeof(*),a::ArrayBlock,b::ArrayBlock)
  c1 = CachedArray(a*b)
  c2 = return_cache(unwrap_cached_array,c1)
  return c1, c2
end

function evaluate!(cache,::typeof(*),a::ArrayBlock,b::ArrayBlock)
  c1, c2 = cache
  setsize_op!(*,c1,a,b)
  c = evaluate!(c2,unwrap_cached_array,c1)
  mul!(c,a,b)
  c
end

# MulAddMap: Cached version of `mul!(d,a,b,α,β)`

struct MulAddMap{T} <: Map
  α::T
  β::T
end

function return_cache(k::MulAddMap,a,b,c)
  d = a*b+c
  CachedArray(d)
end

function evaluate!(cache,k::MulAddMap,a,b,c)
  setsize!(cache,size(c))
  d = cache.array
  copyto!(d,c)
  iszero(k.α) && isone(k.β) && return d
  mul!(d,a,b,k.α,k.β)
  d
end

function return_value(k::MulAddMap,a::ArrayBlock,b::ArrayBlock,c::ArrayBlock)
  x = return_value(*,a,b)
  return_value(+,x,c)
end

function return_cache(k::MulAddMap,a::ArrayBlock,b::ArrayBlock,c::ArrayBlock)
  c1 = CachedArray(a*b+c)
  c2 = return_cache(unwrap_cached_array,c1)
  (c1,c2)
end

function evaluate!(cache,k::MulAddMap,a::ArrayBlock,b::ArrayBlock,c::ArrayBlock)
  c1,c2 = cache
  setsize_op!(copy,c1,c)
  setsize_op!(*,c1,a,b)
  d = evaluate!(c2,unwrap_cached_array,c1)
  copyto!(d,c)
  iszero(k.α) && isone(k.β) && return d
  mul!(d,a,b,k.α,k.β)
  d
end

# Assembly Maps: AddEntriesMap and TouchEntriesMap

struct AddEntriesMap{F} <: Map
  combine::F
end

function evaluate!(cache,k::AddEntriesMap,A,v,i,j)
  add_entries!(k.combine,A,v,i,j)
end

function evaluate!(cache,k::AddEntriesMap,A,v,i)
  add_entries!(k.combine,A,v,i)
end

struct TouchEntriesMap <: Map end

function evaluate!(cache,k::TouchEntriesMap,A,v,i,j)
  add_entries!(+,A,nothing,i,j)
end

function evaluate!(cache,k::TouchEntriesMap,A,v,i)
  add_entries!(+,A,nothing,i)
end

for T in (:AddEntriesMap,:TouchEntriesMap)
  @eval begin

    function return_cache(k::$T,A,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
      qs = findall(v.touched)
      i, j = Tuple(first(qs))
      cij = return_cache(k,A,v.array[i,j],I.array[i],J.array[j])
      ni,nj = size(v.touched)
      cache = Matrix{typeof(cij)}(undef,ni,nj)
      for j in 1:nj
        for i in 1:ni
          if v.touched[i,j]
            cache[i,j] = return_cache(k,A,v.array[i,j],I.array[i],J.array[j])
          end
        end
      end
      cache
    end

    function evaluate!(cache, k::$T,A,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
      ni,nj = size(v.touched)
      for j in 1:nj
        for i in 1:ni
          if v.touched[i,j]
            evaluate!(cache[i,j],k,A,v.array[i,j],I.array[i],J.array[j])
          end
        end
      end
    end

    function return_cache(k::$T,A,v::VectorBlock,I::VectorBlock)
      qs = findall(v.touched)
      i = first(qs)
      ci = return_cache(k,A,v.array[i],I.array[i])
      ni = length(v.touched)
      cache = Vector{typeof(ci)}(undef,ni)
      for i in 1:ni
        if v.touched[i]
          cache[i] = return_cache(k,A,v.array[i],I.array[i])
        end
      end
      cache
    end

    function evaluate!(cache, k::$T,A,v::VectorBlock,I::VectorBlock)
      ni = length(v.touched)
      for i in 1:ni
        if v.touched[i]
          evaluate!(cache[i],k,A,v.array[i],I.array[i])
        end
      end
    end
  end

  for MT in (:MatrixBlock,:MatrixBlockView)
    Aij = (MT == :MatrixBlock) ? :(A.array[i,j]) : :(A[i,j])
    @eval begin

      function return_cache(k::$T,A::$MT,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
        qs = findall(v.touched)
        i, j = Tuple(first(qs))
        cij = return_cache(k,$Aij,v.array[i,j],I.array[i],J.array[j])
        ni,nj = size(v.touched)
        cache = Matrix{typeof(cij)}(undef,ni,nj)
        for j in 1:nj
          for i in 1:ni
            if v.touched[i,j]
              cache[i,j] = return_cache(k,$Aij,v.array[i,j],I.array[i],J.array[j])
            end
          end
        end
        cache
      end

      function evaluate!(cache,k::$T,A::$MT,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
        ni,nj = size(v.touched)
        for j in 1:nj
          for i in 1:ni
            if v.touched[i,j]
              evaluate!(cache[i,j],k,$Aij,v.array[i,j],I.array[i],J.array[j])
            end
          end
        end
      end

    end # @eval
  end # for MT

  for VT in (:VectorBlock,:VectorBlockView)
    Ai = (VT == :VectorBlock) ? :(A.array[i]) : :(A[i])
    @eval begin

      function return_cache(k::$T,A::$VT,v::VectorBlock,I::VectorBlock)
        qs = findall(v.touched)
        i = first(qs)
        ci = return_cache(k,$Ai,v.array[i],I.array[i])
        ni = length(v.touched)
        cache = Vector{typeof(ci)}(undef,ni)
        for i in 1:ni
          if v.touched[i]
            cache[i] = return_cache(k,$Ai,v.array[i],I.array[i])
          end
        end
        cache
      end

      function evaluate!(cache, k::$T,A::$VT,v::VectorBlock,I::VectorBlock)
        ni = length(v.touched)
        for i in 1:ni
          if v.touched[i]
            evaluate!(cache[i],k,$Ai,v.array[i],I.array[i])
          end
        end
      end

    end # @eval
  end # for VT
end # end for T
