
# Target block layout in multi-field computations
#
# u [c][1 , 1][1 ,bu][np,1 ,nu]
# v [c][1 ,  ][bv,  ][np,nv,  ]
# u [f][lc][1 ,lc][1 ,bu][np,1 ,nu]
# v [f][lc][lc,  ][bv,  ][np,nv,  ]
# u⁺[f][1 ,lc][1 ,bu][np,1 ,nu]
# v⁺[f][lc,  ][bv,  ][np,nv,  ]
# u [c][lf][1 ,lF][1 ,bu][np,1 ,nu]
# v [c][lf][lF,  ][bv,  ][np,nv,  ]
#
# λ [f][1 , 1][1 ,bλ][np,1 ,nλ]
# μ [f][1 ,  ][bμ,  ][np,nμ,  ]
# λ [c][lf][1 ,lf][1 ,bλ][np,1 ,nλ]
# μ [c][lf][lf,  ][bμ,  ][np,nμ,  ]
#
# uh [c][np]
# uh [f][lc][np]
# uh⁺[f][np]
# uh [c][lf][np]

"""
    struct ArrayBlock{A,N}
      array::Array{A,N}
      touched::Array{Bool,N}
    end

Block-wise storage of arrays, where `touched` indicates which blocks are active. Accessing 
a non-touched block returns nothing. 

`ArrayBlock` is mostly used for multi-field and skeleton computations and assembly, 
where each block corresponds to a field or plus/minus respectively.
Blocks might be nested in the case of multi-field computations on skeletons.

"""
struct ArrayBlock{A,N}
  array::Array{A,N}
  touched::Array{Bool,N}
  function ArrayBlock(array::Array{A,N},touched::Array{Bool,N}) where {A,N}
    @check size(array) == size(touched)
    new{A,N}(array,touched)
  end
end

get_array(b::ArrayBlock) = b.array

const VectorBlock = ArrayBlock{A,1} where A
const MatrixBlock = ArrayBlock{A,2} where A

Base.axes(b::ArrayBlock,i) = axes(b.array,i)
Base.size(b::ArrayBlock) = size(b.array)
Base.size(b::ArrayBlock,i) = size(b.array,i)
Base.length(b::ArrayBlock) = length(b.array)
Base.eltype(::Type{<:ArrayBlock{A}}) where A = A
Base.eltype(b::ArrayBlock{A}) where A = A
Base.ndims(b::ArrayBlock{A,N}) where {A,N} = N
Base.ndims(::Type{ArrayBlock{A,N}}) where {A,N} = N
function Base.getindex(b::ArrayBlock,i...)
  if !b.touched[i...]
    return nothing
  end
  b.array[i...]
end
function Base.setindex!(b::ArrayBlock,v,i...)
 @check b.touched[i...] "Only touched entries can be set"
 b.array[i...] = v
end
function Base.show(io::IO,o::ArrayBlock)
  print(io,"ArrayBlock($(o.array), $(o.touched))")
end
function Base.show(io::IO,k::MIME"text/plain",o::ArrayBlock)
  function _sub_block_str(a)
    function _nice_typeof(v)
      "$v"
    end
    function _nice_typeof(::Type{<:VectorBlock{A}}) where A
      "VectorBlock{$(_nice_typeof(A))}"
    end
    function _nice_typeof(::Type{<:MatrixBlock{A}}) where A
      "MatrixBlock{$(_nice_typeof(A))}"
    end
    if isa(a,AbstractArray) || isa(a,ArrayBlock)
      n = replace(_nice_typeof(typeof(a)),"Gridap.Fields."=>"")
      if ndims(a) == 1
        "$(length(a))-element $n"
      else
        s = replace(replace(replace("$(size(a))",", "=>"x"),"("=>""),")"=>"")
        "$s $n"
      end
    else
      "$a"
    end
  end
  println(io,"$(_sub_block_str(o)):")
  lis = LinearIndices(o.array)
  for i in CartesianIndices(o.array)
    if ndims(o) == 1
      s = "[$(lis[i])]"
    else
      s = replace(replace("$(Tuple(i))","("=>"["),")"=>"]")
    end
    if lis[i] != 1
      print(io,"\n")
    end
    print(io," $s = $(_sub_block_str(o[i]))")
  end
end

function Base.similar(x::ArrayBlock{A,N}, T) where {A,N}
  touched = copy(x.touched)
  B = typeof(similar(testvalue(A),T))
  array = Array{B,N}(undef, size(touched))
  for I in eachindex(touched)
    if touched[I]
      array[I] = similar(x.array[I],T)
    end
  end
  ArrayBlock(array,touched)
end

function Base.:≈(a::AbstractArray{<:ArrayBlock},b::AbstractArray{<:ArrayBlock})
  all(z->z[1]≈z[2],zip(a,b))
end

function Base.:≈(a::ArrayBlock,b::ArrayBlock)
  if size(a) != size(b)
    return false
  end
  if a.touched != b.touched
    return false
  end
  for i in eachindex(a.array)
    if a.touched[i]
      if !(a.array[i] ≈ b.array[i])
        return false
      end
    end
  end
  true
end

function Base.:(==)(a::ArrayBlock,b::ArrayBlock)
  if size(a) != size(b)
    return false
  end
  for i in eachindex(a.array)
    if a.touched[i] && b.touched[i]
      if a.array[i] != b.array[i]
        return false
      end
    elseif a.touched[i]
      return false
    elseif b.touched[i]
      return false
    end
  end
  true
end

Base.copy(a::ArrayBlock) = ArrayBlock(copy(a.array),copy(a.touched))
Base.eachindex(a::ArrayBlock) = eachindex(a.array)

entry_type(a) = eltype(a)
entry_type(a::ArrayBlock) = entry_type(typeof(a))
entry_type(::Type{ArrayBlock{A,N}}) where {A,N} = entry_type(A)

fill_entries!(a,value) = fill!(a,value)
function fill_entries!(a::ArrayBlock,value)
  for i in eachindex(a.array)
    if a.touched[i]
      fill_entries!(a.array[i],value)
    end
  end
  return a
end

function Base.copyto!(d::ArrayBlock,c::ArrayBlock)
  z = zero(entry_type(d))
  for i in eachindex(c.array)
    if c.touched[i]
      copyto!(d.array[i],c.array[i])
    elseif d.touched[i]
      fill_entries!(d.array[i],z)
    end
  end
  d
end

function testitem(f::ArrayBlock{A}) where A
  @notimplementedif !isconcretetype(A)
  i = findall(f.touched)
  if length(i) != 0
    f.array[i[1]]
  else
    testvalue(A)
  end
end

function testvalue(::Type{ArrayBlock{A,N}}) where {A,N}
  s = ntuple(i->0,Val(N))
  array = Array{A,N}(undef,s)
  touched = Array{Bool,N}(undef,s)
  ArrayBlock(array,touched)
end

# CachedArray methods

function CachedArray(a::ArrayBlock)
  ai = testitem(a)
  ci = CachedArray(ai)
  array = Array{typeof(ci),ndims(a)}(undef,size(a))
  for i in eachindex(a.array)
    if a.touched[i]
      array[i] = CachedArray(a.array[i])
    end
  end
  ArrayBlock(array,a.touched)
end

function setsize_op!(::typeof(copy),a::ArrayBlock,b::ArrayBlock)
  @check size(a) == size(b)
  for i in eachindex(b.array)
    if b.touched[i]
      setsize_op!(copy,a.array[i],b.array[i])
    end
  end
  return a
end

function setsize_op!(::typeof(*),c::VectorBlock,a::MatrixBlock,b::VectorBlock)
  ni, nj = size(a)
  @check length(b) == nj
  @check length(c) == ni
  for i in 1:ni
    for j in 1:nj
      if a.touched[i,j] && b.touched[j]
        setsize_op!(*,c.array[i],a.array[i,j],b.array[j])
      end
    end
  end
  c
end

function setsize_op!(::typeof(*),c::MatrixBlock,a::MatrixBlock,b::MatrixBlock)
  ni,nk = size(a)
  nk2,nj = size(b)
  @check nk == nk2
  @check size(c) == (ni,nj)
  for i in 1:ni
    for j in 1:nj
      for k in 1:nk
        if a.touched[i,k] && b.touched[k,j]
          setsize_op!(*,c.array[i,j],a.array[i,k],b.array[k,j])
        end
      end
    end
  end
  c
end

function unwrap_cached_array(a::CachedArray)
  a.array
end

function unwrap_cached_array(a::ArrayBlock)
  cache = return_cache(unwrap_cached_array,a)
  evaluate!(cache,unwrap_cached_array,a)
end

function return_cache(::typeof(unwrap_cached_array),a::ArrayBlock)
  ai = testitem(a)
  ci = return_cache(unwrap_cached_array,ai)
  ri = evaluate!(ci,unwrap_cached_array,ai)
  c = Array{typeof(ci),ndims(a)}(undef,size(a))
  array = Array{typeof(ri),ndims(a)}(undef,size(a))
  for i in eachindex(a.array)
    if a.touched[i]
      c[i] = return_cache(unwrap_cached_array,a.array[i])
    end
  end
  ArrayBlock(array,a.touched), c
end

function evaluate!(cache,::typeof(unwrap_cached_array),a::ArrayBlock)
  r,c = cache
  for i in eachindex(a.array)
    if a.touched[i]
      r.array[i] = evaluate!(c[i],unwrap_cached_array,a.array[i])
    end
  end
  r
end

# BlockMap

"""
    struct BlockMap{N} <: Map
        size::NTuple{N,Int}
        indices::Vector{CartesianIndex{N}}
    end

A `BlockMap` maps `M = length(indices)` arrays to a single array of N-dimensinal blocks, 
where only the blocks indexed by `indices` are touched and contain the corresponding 
entry of the input arrays.

## Constructors:

    BlockMap(l::Integer,i::Integer) ≡ BlockMap((l,),[CartesianIndex((i,))])
    BlockMap(l::Integer,inds::Vector{<:Integer}) ≡ BlockMap((l,),[CartesianIndex((i,)) for i in inds])
    BlockMap(s::NTuple,i::Integer) ≡ BlockMap(s,[CartesianIndex(i)])
    BlockMap(s::NTuple,inds::Vector{<:NTuple}) ≡ BlockMap(s,[CartesianIndex(i) for i in inds])

## Usage:

```julia
  lazy_map(BlockMap(2,[1,2]),a,b)
```
"""
struct BlockMap{N} <: Map
  size::NTuple{N,Int}
  indices::Vector{CartesianIndex{N}}
end

function BlockMap(l::Integer,i::Integer)
  s = (l,)
  ci = CartesianIndex((i,))
  BlockMap(s,[ci])
end

function BlockMap(l::Integer,inds::Vector{<:Integer})
  s = (l,)
  cis = [CartesianIndex((i,)) for i in inds]
  BlockMap(s,cis)
end

function BlockMap(s::NTuple,i::Integer)
  cis = CartesianIndices(s)
  ci = cis[i]
  BlockMap(s,[ci])
end

function BlockMap(s::NTuple,inds::Vector{<:NTuple})
  cis = [CartesianIndex(i) for i in inds]
  BlockMap(s,cis)
end

return_type(k::BlockMap{N},a::A...) where {A,N} = ArrayBlock{A,N}

function return_cache(k::BlockMap{N},a::A...) where {A,N}
  array = Array{A,N}(undef,k.size)
  touched = fill(false,k.size)
  for (t,i) in enumerate(k.indices)
    array[i] = a[t]
    touched[i] = true
  end
  ArrayBlock(array,touched)
end

function evaluate!(cache,k::BlockMap{N},a::A...) where {A,N}
  @check size(cache) == k.size
  for (t,i) in enumerate(k.indices)
    cache.array[i] = a[t]
  end
  cache
end

# ZeroBlockMap

struct ZeroBlockMap <: Map end

function return_cache(::ZeroBlockMap,a::AbstractArray,b::AbstractArray)
  CachedArray(similar(a,eltype(a),size(b)))
end

function evaluate!(cache,::ZeroBlockMap,a,b::AbstractArray)
  setsize!(cache,size(b))
  r = cache.array
  fill!(r,zero(eltype(r)))
  r
end

function return_cache(::ZeroBlockMap,a::ArrayBlock,b::ArrayBlock)
  A = eltype(a)
  N = ndims(b)
  array = Array{A,N}(undef,size(b))
  touched = fill(false,size(b))
  ArrayBlock(array,touched)
end

function evaluate!(cache,::ZeroBlockMap,a,b::ArrayBlock)
  @check size(cache) == size(b)
  cache
end

"""
    struct MergeBlockMap{N,M} <: Map
        size::NTuple{N,Int}
        indices::Vector{Vector{Tuple{CartesianIndex{M},CartesianIndex{N}}}}
    end

A `MergeBlockMap` create a single array of N-dimensional blocks from `L=length(indices)` 
input arrays of M-dimensional blocks. 

For the l-th input array `a_l`, the vector of tuples in `indices[l]` contains the 
mapping between the indices of the blocks in `a_l` and the indices of the blocks in
the output array, i.e `a_out[P[2]] = a_l[P[1]] ∀ P ∈ indices[l], ∀ l `.
"""
struct MergeBlockMap{N,M} <: Map
  size::NTuple{N,Int}
  indices::Vector{Vector{Tuple{CartesianIndex{M},CartesianIndex{N}}}}
end

function return_cache(k::MergeBlockMap{N,M},a::ArrayBlock{A,M}...) where {A,N,M}
  array = Array{A,N}(undef,k.size)
  touched = fill(false,k.size)
  for (t,I) in enumerate(k.indices)
    at = a[t]
    for (j,i) in I
      if at.touched[j]
        array[i] = at[j]
        touched[i] = true
      end
    end
  end
  ArrayBlock(array,touched)
end

function evaluate!(cache,k::MergeBlockMap{N,M},a::ArrayBlock{A,M}...) where {A,N,M}
  @check size(cache) == k.size
  for (t,I) in enumerate(k.indices)
    at = a[t]
    for (j,i) in I
      if at.touched[j]
        cache.array[i] = at[j]
      end
    end
  end
  cache
end

"""
    struct BlockBroadcasting{F} <: Map
      f::F
    end

A `BlockBroadcasting` broadcasts a map `f` block-wise over the input arrays.
"""
struct BlockBroadcasting{F} <: Map
  f::F
end

function return_value(k::BlockBroadcasting,a::ArrayBlock{A,N},b::ArrayBlock...) where {A,N}
  @check all(a.touched == bi.touched for bi in b)

  i = findfirst(a.touched)
  ai = (a.array[i],(bi.array[i] for bi in b)...)
  ri = return_value(k.f,ai...)

  array = Array{typeof(ri),N}(undef,size(a))
  for i in eachindex(a.array)
    if a.touched[i]
      ai = (a.array[i],(bi.array[i] for bi in b)...)
      array[i] = return_value(k.f,ai...)
    end
  end

  return ArrayBlock(array,a.touched)
end

function return_cache(k::BlockBroadcasting,a::ArrayBlock{A,N},b::ArrayBlock...) where {A,N}
  @check all(a.touched == bi.touched for bi in b)

  i = findfirst(a.touched)
  ai = (a.array[i],(bi.array[i] for bi in b)...)
  ci = return_cache(k.f,ai...)
  ri = evaluate!(ci,k.f,ai...)

  array = Array{typeof(ri),N}(undef,size(a))
  caches = Array{typeof(ci),N}(undef,size(a))
  for i in eachindex(a.array)
    if a.touched[i]
      ai = (a.array[i],(bi.array[i] for bi in b)...)
      caches[i] = return_cache(k.f,ai...)
    end
  end

  return ArrayBlock(array,a.touched), caches
end

function evaluate!(cache,k::BlockBroadcasting,a::ArrayBlock{A,N},b::ArrayBlock...) where {A,N}
  r, c = cache
  @check r.touched == a.touched
  @check all(a.touched == bi.touched for bi in b)
  
  for i in eachindex(a.array)
    if a.touched[i]
      ai = (a.array[i],(bi.array[i] for bi in b)...)
      r.array[i] = evaluate!(c[i],k.f,ai...)
    end
  end

  return r
end

# Multiplication on ArrayBlocks

function Base.:*(a::Number,b::ArrayBlock)
  bi = testitem(b)
  ci = a*bi
  array = Array{typeof(ci),ndims(b)}(undef,size(b))
  for i in eachindex(b.array)
    if b.touched[i]
      array[i] = a*b.array[i]
    end
  end
  ArrayBlock(array,b.touched)
end

Base.:*(a::ArrayBlock,b::Number) = b*a

function Base.:*(a::ArrayBlock{A,2},b::ArrayBlock{B,1}) where {A,B}
  @check size(a.array,2) == size(b.array,1)
  C = typeof(testvalue(A)*testvalue(B))
  array = Vector{C}(undef,size(a.array,1))
  touched = fill(false,size(a.array,1))
  ni,nj = size(a.array)
  for i in 1:ni
    for j in 1:nj
      if a.touched[i,j] && b.touched[j]
        if !touched[i]
          array[i] = a.array[i,j]*b.array[j]
          touched[i] = true
        elseif isa(C,AbstractArray)
          mul!(array[i],a.array[i,j],b.array[j],1,1)
        else
          array[i] += a.array[i,j]*b.array[j]
        end
      end
    end
  end
  ArrayBlock(array,touched)
end

function Base.:*(a::ArrayBlock{A,1},b::ArrayBlock{B,2}) where {A,B}
  @check size(b.array,1) == 1
  ai = testvalue(A)
  bi = testvalue(B)
  ri = ai*bi
  ni, nj = size(a.array,1), size(b.array,2)
  array = Matrix{typeof(ri)}(undef,ni,nj)
  touched = fill(false,ni,nj)
  for i in 1:ni
    for j in 1:nj
      if a.touched[i] && b.touched[1,j]
        array[i,j]   = a.array[i]*b.array[1,j]
        touched[i,j] = true
      end
    end
  end
  ArrayBlock(array,touched)
end

function Base.:*(a::ArrayBlock{A,2},b::ArrayBlock{B,2}) where {A,B}
  @check size(a.array,2) == size(b.array,1)
  C = typeof(testvalue(A)*testvalue(B))
  ni = size(a.array,1)
  nj = size(b.array,2)
  nk = size(a.array,2)
  array = Matrix{C}(undef,(ni,nj))
  touched = fill(false,(ni,nj))
  ni,nj = size(a.array)
  for i in 1:ni
    for j in 1:nj
      for k in 1:nk
        if a.touched[i,k] && b.touched[k,j]
          if !touched[i,j]
            array[i,j] = a.array[i,k]*b.array[k,j]
            touched[i,j] = true
          elseif isa(C,AbstractArray)
            mul!(array[i,j],a.array[i,k],b.array[k,j],1,1)
          else
            array[i,j] += a.array[i,k]*b.array[k,j]
          end
        end
      end
    end
  end
  ArrayBlock(array,touched)
end

# LinearAlgebra-related methods

function LinearAlgebra.rmul!(a::ArrayBlock,β)
  for i in eachindex(a.touched)
    if a.touched[i]
      rmul!(a.array[i],β)
    end
  end
end

function LinearAlgebra.mul!(
  c::ArrayBlock{C,1} where C,
  a::ArrayBlock{A,2} where A,
  b::ArrayBlock{B,1} where B,
  α::Number,β::Number)

  ni, nj = size(a)
  @check nj == size(b.array,1)
  for i in 1:ni
    if β!=1 && c.touched[i]
      rmul!(c.array[i],β)
    end
    for j in 1:nj
      if a.touched[i,j] && b.touched[j]
        ci = c.array[i]
        aij = a.array[i,j]
        bj = b.array[j]
        mul!(ci,aij,bj,α,1)
      end
    end
  end
  c
end

function LinearAlgebra.mul!(
  c::ArrayBlock{C,2} where C,
  a::ArrayBlock{A,2} where A,
  b::ArrayBlock{B,2} where B,
  α::Number,β::Number)

  ni, nk = size(a.array)
  nj = size(b.array,2)
  @check nk == size(b.array,1)
  @check (ni,nj) == size(c.array)

  for i in 1:ni
    for j in 1:nj
      if β!=1 && c.touched[i,j]
        rmul!(c.array[i,j],β)
      end
      for k in 1:nk
        if a.touched[i,k] && b.touched[k,j]
          cij = c.array[i,j]
          aik = a.array[i,k]
          bkj = b.array[k,j]
          mul!(cij,aik,bkj,α,1)
        end
      end
    end
  end
  c
end

# Broadcasting on ArrayBlocks

function return_value(k::Broadcasting,f::ArrayBlock{A,N}) where {A,N}
  fi = testitem(f)
  fix = return_value(k,fi)
  g = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      g[i] = return_value(k,f.array[i])
    end
  end
  ArrayBlock(g,f.touched)
end

function return_cache(k::Broadcasting,f::ArrayBlock{A,N}) where {A,N}
  fi = testitem(f)
  li = return_cache(k,fi)
  fix = evaluate!(li,k,fi)
  l = Array{typeof(li),N}(undef,size(f.array))
  g = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      l[i] = return_cache(k,f.array[i])
    end
  end
  ArrayBlock(g,f.touched),l
end

function evaluate!(cache,k::Broadcasting,f::ArrayBlock)
  g,l = cache
  @check g.touched == f.touched
  for i in eachindex(f.array)
    if f.touched[i]
      g.array[i] = evaluate!(l[i],k,f.array[i])
    end
  end
  g
end

function return_value(k::Broadcasting{typeof(*)},f::Number,g::ArrayBlock)
  gi = testitem(g)
  hi = return_value(k,f,gi)
  array = Array{typeof(hi),ndims(g.array)}(undef,size(g.array))
  for i in eachindex(g.array)
    if g.touched[i]
      array[i] = return_value(k,f,g.array[i])
    end
  end
  ArrayBlock(array,g.touched)
end

function return_cache(k::Broadcasting{typeof(*)},f::Number,g::ArrayBlock)
  gi = testitem(g)
  ci = return_cache(k,f,gi)
  hi = evaluate!(ci,k,f,gi)
  array = Array{typeof(hi),ndims(g.array)}(undef,size(g.array))
  c = Array{typeof(ci),ndims(g.array)}(undef,size(g.array))
  for i in eachindex(g.array)
    if g.touched[i]
      c[i] = return_cache(k,f,g.array[i])
    end
  end
  ArrayBlock(array,g.touched), c
end

function evaluate!(cache,k::Broadcasting{typeof(*)},f::Number,g::ArrayBlock)
  r, c = cache
  for i in eachindex(g.array)
    if g.touched[i]
      r.array[i] = evaluate!(c[i],k,f,g.array[i])
    end
  end
  r
end

function return_value(k::Broadcasting{typeof(*)},f::ArrayBlock,g::Number)
  return_value(k,g,f)
end

function return_cache(k::Broadcasting{typeof(*)},f::ArrayBlock,g::Number)
  return_cache(k,g,f)
end

function evaluate!(cache,k::Broadcasting{typeof(*)},f::ArrayBlock,g::Number)
  evaluate!(cache,k,g,f)
end

# ArrayBlock views

struct ArrayBlockView{A,N,M}
  array::ArrayBlock{A,M}
  block_map::Array{CartesianIndex{M},N}
end

Base.view(a::ArrayBlock{A,M},b::Array{CartesianIndex{M},N}) where {A,M,N} = ArrayBlockView(a,b)
const MatrixBlockView{A} = ArrayBlockView{A,2,2} where A
const VectorBlockView{A} = ArrayBlockView{A,1,1} where A

Base.axes(a::ArrayBlockView,i) = axes(a.block_map,i)
Base.size(a::ArrayBlockView) = size(a.block_map)
Base.length(b::ArrayBlockView) = length(b.block_map)
Base.eltype(::Type{<:ArrayBlockView{A}}) where A = A
Base.eltype(::ArrayBlockView{A}) where A = A
Base.ndims(::ArrayBlockView{A,N}) where {A,N} = N
Base.ndims(::Type{ArrayBlockView{A,N}}) where {A,N} = N
Base.getindex(b::ArrayBlockView,i...) = getindex(b.array,b.block_map[i...])
Base.setindex!(b::ArrayBlockView,v,i...) = setindex!(b.array,v,b.block_map[i...])

Base.copy(a::ArrayBlockView) = ArrayBlockView(copy(a.array),copy(a.block_map))
Base.eachindex(a::ArrayBlockView) = eachindex(a.block_map)

function Base.show(io::IO,o::ArrayBlockView)
  print(io,"ArrayBlockView($(o.array), $(o.block_map))")
end

LinearAlgebra.diag(a::MatrixBlockView) = view(a.array.array, diag(a.block_map))
LinearAlgebra.diag(a::MatrixBlock) = view(a.array,diag(CartesianIndices(a.array)))

entry_type(a::ArrayBlockView) = entry_type(a.array)
fill_entries!(a::ArrayBlockView,value) = fill_entries!(a.array,value)

function Arrays.CachedArray(a::ArrayBlockView)
  ArrayBlockView(CachedArray(a.array),a.block_map)
end

function unwrap_cached_array(a::ArrayBlockView)
  cache = return_cache(unwrap_cached_array,a)
  evaluate!(cache,unwrap_cached_array,a)
end

function return_cache(::typeof(unwrap_cached_array),a::ArrayBlockView)
  cache = return_cache(unwrap_cached_array,a.array)
  array = evaluate!(cache,unwrap_cached_array,a.array)
  return ArrayBlockView(array,a.block_map), cache
end

function evaluate!(cache,::typeof(unwrap_cached_array),a::ArrayBlockView)
  r, c = cache
  evaluate!(c,unwrap_cached_array,a.array)
  return r
end
