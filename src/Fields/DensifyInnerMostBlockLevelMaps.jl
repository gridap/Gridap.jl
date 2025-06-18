"""
    struct DensifyInnerMostBlockLevelMap <: Map
"""
struct DensifyInnerMostBlockLevelMap <: Map
end

function return_cache(k::DensifyInnerMostBlockLevelMap,
                      a::VectorBlock{<:Vector{T}}) where {T}
  @check all(a.touched)
  # Precompute the size of each row block.
  # We are assuming that the size of each row
  # block is going to be equivalent for all cells
  s=0
  brs=Vector{Int}(undef,size(a)[1])
  for i=1:size(a)[1]
    s=s+length(a.array[i])
    brs[i]=length(a.array[i])
  end
  CachedArray(Vector{T}(undef,(s,))),brs
end

function return_cache(k   :: DensifyInnerMostBlockLevelMap,
                      brs :: Vector{<:Integer},
                      a   :: VectorBlock{<:Vector{T}}) where {T}
  @check length(brs)==size(a)[1]
  CachedArray(Vector{T}(undef,(sum(brs),)))
end


function return_cache(k::DensifyInnerMostBlockLevelMap,
                      a::MatrixBlock{<:Vector{T}}) where {T}
  @check _check_preconditions(k,a)
  s=Vector{Int}(undef,2)
  s[1]=0
  s[2]=size(a)[2]
  # Precompute the size of each row block.
  # We are assuming that the size of each row
  # block is going to be equivalent for all cells
  brs=Vector{Int}(undef,size(a)[1])
  for i=1:size(a)[1]
    for j=1:size(a)[2]
      if (a.touched[i,j])
        s[1]=s[1]+length(a.array[i,j])
        brs[i]=length(a.array[i,j])
        break
      end
    end
  end
  CachedArray(Array{T,2}(undef,Tuple(s))),brs
end

function return_cache(k   :: DensifyInnerMostBlockLevelMap,
                      brs :: Vector{<:Integer},
                      a   :: MatrixBlock{<:Vector{T}}) where {T}
  @check length(brs)==size(a)[1]
  CachedArray(Array{T,2}(undef,(sum(brs),size(a)[2])))
end

function _check_preconditions(k::DensifyInnerMostBlockLevelMap,
                              a::MatrixBlock{<:Vector{T}}) where {T}
  # Double check that, for each row, there is at least one non-zero block
  found=false
  for i=1:size(a)[1]
    found=false
    for j=1:size(a)[2]
      if (a.touched[i,j])
        found = true
        break
      end
    end
    if !found
      break
    end
  end
  found
end

function return_cache(k::DensifyInnerMostBlockLevelMap,
                      a::VectorBlock{<:Matrix{T}}) where {T}
  @check all(a.touched)
  s=Vector{Int}(undef,2)
  s[1]=0
  s[2]=size(a.array[1])[2]
  for i=1:size(a)[1]
    s[1]=s[1]+size(a.array[i])[1]
  end
  CachedArray(Array{T,2}(undef,Tuple(s)))
end

function return_cache(k  ::DensifyInnerMostBlockLevelMap,
                      brs::Vector{<:Integer},
                      cs ::Integer,
                      a  ::VectorBlock{<:Matrix{T}}) where {T}
  CachedArray(Array{T,2}(undef,(sum(brs),cs)))
end


function return_cache(k::DensifyInnerMostBlockLevelMap,
                      a::MatrixBlock{<:Matrix{T}}) where {T}
  @check _check_preconditions(k,a)
  s=Vector{Int}(undef,2)
  s[1]=0
  s[2]=0
  # Precompute the size of each row/col block.
  # We are assuming that the size of each row/col
  # block is going to be equivalent for all cells
  brs=Vector{Int}(undef,size(a)[1])
  bcs=Vector{Int}(undef,size(a)[2])
  # Row traversal
  for i=1:size(a)[1]
    for j=1:size(a)[2]
      if (a.touched[i,j])
        s[1]=s[1]+size(a.array[i,j])[1]
        brs[i]=size(a.array[i,j])[1]
        break
      end
    end
  end
  #Column traversal
  for j=1:size(a)[2]
     for i=1:size(a)[1]
       if (a.touched[i,j])
        s[2]=s[2]+size(a.array[i,j])[2]
        bcs[j]=size(a.array[i,j])[2]
        break
       end
     end
  end
  CachedArray(Array{T,2}(undef,Tuple(s))),brs,bcs
end

function _check_preconditions(k::DensifyInnerMostBlockLevelMap,
                              a::MatrixBlock{<:Matrix{T}}) where {T}
  # Double check that, for each row,
  # there is at least one non-zero block
  found=false
  for i=1:size(a)[1]
    found=false
    for j=1:size(a)[2]
      if (a.touched[i,j])
        found = true
        break
      end
    end
    if !found
      break
    end
  end
  # Double check that, for each col,
  # there is at least one non-zero block
  if (found)
    for j=1:size(a)[2]
      found=false
      for i=1:size(a)[1]
        if (a.touched[i,j])
          found = true
          break
        end
      end
      if !found
        break
      end
    end
  end
  found
end

function return_cache(k::DensifyInnerMostBlockLevelMap,
                      brs::Vector{<:Integer},
                      bcs::Vector{<:Integer},
                      a::MatrixBlock{<:Matrix{T}}) where {T}
  CachedArray(Array{T,2}(undef,(sum(brs),sum(bcs))))
end


function evaluate!(cache,
                   k::DensifyInnerMostBlockLevelMap,
                   a::VectorBlock{<:Vector{T}}) where {T}
  cache_array, brs = cache
  _evaluate!(cache_array,brs,a)
end

function evaluate!(cache,
                   k   :: DensifyInnerMostBlockLevelMap,
                   brs :: Vector{<:Integer},
                   a   :: VectorBlock{<:Vector{T}}) where {T}
  _evaluate!(cache,brs,a)
end

function _evaluate!(cache,
                    brs::Vector{<:Integer},
                    a::VectorBlock{<:Vector{T}}) where {T}
  output  = cache.array
  output .= zero(eltype(output))
  current_i=1
  for i=1:size(a)[1]
    range = current_i:current_i+brs[i]-1
    if (a.touched[i])
      output[range] = a.array[i]
    end
    current_i = current_i + length(range)
  end
  output
end

function evaluate!(cache,
                   k::DensifyInnerMostBlockLevelMap,
                   a::MatrixBlock{<:Vector{T}}) where {T}
  cache_array, brs = cache
  _evaluate!(cache_array,brs,a)
end

function evaluate!(cache,
                   k   :: DensifyInnerMostBlockLevelMap,
                   brs :: Vector{<:Integer},
                   a   :: MatrixBlock{<:Vector{T}}) where {T}
  _evaluate!(cache,brs,a)
end

function _evaluate!(cache,
                    brs::Vector{<:Integer},
                    a::MatrixBlock{<:Vector{T}}) where {T}
  output  = cache.array
  output .= zero(eltype(output))
  for j=1:size(a)[2]
     current_i=1
     for i=1:size(a)[1]
        range = current_i:current_i+brs[i]-1
        if (a.touched[i,j])
          output[range,j] = a.array[i,j]
        end
        current_i = current_i + length(range)
     end
  end
  output
end

function evaluate!(cache,
                   k::DensifyInnerMostBlockLevelMap,
                   a::MatrixBlock{<:Matrix{T}}) where {T}
  cache_array, brs, bcs = cache
  _evaluate!(cache_array, brs, bcs,a)
end

function evaluate!(cache,
  k::DensifyInnerMostBlockLevelMap,
  brs::Vector{<:Integer},
  bcs::Vector{<:Integer},
  a::MatrixBlock{<:Matrix{T}}) where {T}
  _evaluate!(cache, brs, bcs, a)
end

function _evaluate!(cache_array, brs, bcs, a)
  output  = cache_array.array
  output .= zero(eltype(output))
  current_j=1
  for j=1:size(a)[2]
    current_i=1
    range_j = current_j:current_j+bcs[j]-1
    for i=1:size(a)[1]
      range_i = current_i:current_i+brs[i]-1
      if (a.touched[i,j])
         output[range_i,range_j] = a.array[i,j]
      end
      current_i = current_i + brs[i]
    end
    current_j = current_j + bcs[j]
  end
  output
end

function evaluate!(cache,
                   k::DensifyInnerMostBlockLevelMap,
                   a::VectorBlock{<:Matrix{T}}) where {T}
  @check all(a.touched)
  output = cache.array
  current_i=1
  n=size(a.array[1])[2]
  for i=1:size(a)[1]
    range=current_i:current_i+size(a.array[i])[1]-1
    output[range,1:n] = a.array[i]
    current_i=current_i+length(range)
  end
  output
end

function return_cache(k::DensifyInnerMostBlockLevelMap,
      a::ArrayBlock{<:ArrayBlock{T,M} where {T,M},N}) where {N}
    cache_touched=a.touched
    i=findfirst(isone, cache_touched)
    cache_block=return_cache(k,a.array[i])
    cache_array=Array{typeof(cache_block),N}(undef,size(a))
    output_array=Array{return_type(k,a.array[i]),N}(undef,size(a))
    linds=LinearIndices(size(a))
    cinds=CartesianIndices(size(a))
    while i != nothing
      cache_array[i]=cache_block
      if (linds[i]+1 <= length(linds))
        i=findnext(isone, cache_touched, cinds[linds[i]+1])
        if (i!=nothing)
          cache_block=return_cache(k,a.array[i])
        end
      else
        i=nothing
      end
    end
    ArrayBlock(cache_array,cache_touched),
    ArrayBlock(output_array,cache_touched),
    linds,
    cinds
end

function evaluate!(cache,
    k::DensifyInnerMostBlockLevelMap,
    a::ArrayBlock{<:ArrayBlock{T,M} where {T,M},N}) where {N}
    cache_array, output_array, linds, cinds = cache
    @check cache_array.touched == a.touched
    i=findfirst(isone, cache_array.touched)
    while i != nothing
      output_array.array[i]=evaluate!(cache_array.array[i],k,a.array[i])
      if (linds[i]+1 <= length(linds))
        i=findnext(isone, cache_array.touched, cinds[linds[i]+1])
      else
        i=nothing
      end
    end
    output_array
end
