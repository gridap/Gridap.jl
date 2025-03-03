
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


struct ArrayBlock{A,N}
  array::Array{A,N}
  touched::Array{Bool,N}
  function ArrayBlock(array::Array{A,N},touched::Array{Bool,N}) where {A,N}
    @check size(array) == size(touched)
    new{A,N}(array,touched)
  end
end

const VectorBlock = ArrayBlock{A,1} where A
const MatrixBlock = ArrayBlock{A,2} where A

Base.axes(b::ArrayBlock,i) = axes(b.array,i)
Base.size(b::ArrayBlock) = size(b.array)
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

<<<<<<< HEAD
Arrays.get_array(b::ArrayBlock) = b.array
=======
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
>>>>>>> 6db1059ae004afaa5ba9764ed865e27872bfc5ac

function Arrays.testitem(f::ArrayBlock{A}) where A
  @notimplementedif !isconcretetype(A)
  i = findall(f.touched)
  if length(i) != 0
    f.array[i[1]]
  else
    testvalue(A)
  end
end

function Arrays.testvalue(::Type{ArrayBlock{A,N}}) where {A,N}
  s = ntuple(i->0,Val(N))
  array = Array{A,N}(undef,s)
  touched = Array{Bool,N}(undef,s)
  ArrayBlock(array,touched)
end

function Arrays.CachedArray(a::ArrayBlock)
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

function lazy_map(k::BlockMap,a::LazyArray{<:Fill{<:Broadcasting{typeof(∘)}}})
  b = a.args[1]
  ϕ = a.args[2]
  c = lazy_map(k,b)
  lazy_map(Broadcasting(∘),c,ϕ)
end

function lazy_map(k::BlockMap,a::MemoArray)
  b = a.parent
  c = lazy_map(k,b)
  MemoArray(c)
end

function lazy_map(k::BlockMap,a::LazyArray{<:Fill{typeof(transpose)}})
  @notimplementedif k.size[1] != 1
  lis = LinearIndices(k.size)
  k2 = BlockMap(k.size[2],[ lis[i] for i in k.indices])
  b = a.args[1]
  c = lazy_map(k2,b)
  lazy_map(transpose,c)
end

function lazy_map(::typeof(evaluate),a::LazyArray{<:Fill{<:BlockMap}},x::AbstractArray)
  args = map(i->lazy_map(evaluate,i,x),a.args)
  k = a.maps.value
  lazy_map(k,args...)
end

# This lazy_map is triggered from function gradient(a::CellField) with optimization
# purposes. See https://github.com/gridap/Gridap.jl/pull/638 for more details.
function lazy_map(k::Broadcasting{typeof(gradient)},a::LazyArray{<:Fill{<:BlockMap}})
  args = map(i->lazy_map(k,i),a.args)
  bm = a.maps.value
  lazy_map(bm,args...)
end

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

struct BlockBroadcasting{F} <: Map
  f::F
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

function return_cache(f::ArrayBlock{A,N},x) where {A,N}
  fi = testitem(f)
  li = return_cache(fi,x)
  fix = evaluate!(li,fi,x)
  l = Array{typeof(li),N}(undef,size(f.array))
  g = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      l[i] = return_cache(f.array[i],x)
    end
  end
  ArrayBlock(g,f.touched),l
end

function evaluate!(cache,f::ArrayBlock,x)
  g,l = cache
  @check g.touched == f.touched
  for i in eachindex(f.array)
    if f.touched[i]
      g.array[i] = evaluate!(l[i],f.array[i],x)
    end
  end
  g
end

function linear_combination(u::ArrayBlock,f::ArrayBlock)
  i::Int = findfirst(f.touched)
  fi = f.array[i]
  ui = u.array[i]
  ufi = linear_combination(ui,fi)
  g = Vector{typeof(ufi)}(undef,length(f.touched))
  for i in eachindex(f.touched)
    if f.touched[i]
      g[i] = linear_combination(u.array[i],f.array[i])
    end
  end
  ArrayBlock(g,f.touched)
end

function return_cache(k::LinearCombinationMap,u::ArrayBlock,fx::ArrayBlock)
  i::Int = findfirst(fx.touched)
  fxi = fx.array[i]
  ui = u.array[i]
  li = return_cache(k,ui,fxi)
  ufxi = evaluate!(li,k,ui,fxi)
  l = Vector{typeof(li)}(undef,size(fx.array))
  g = Vector{typeof(ufxi)}(undef,size(fx.array))
  for i in eachindex(fx.array)
    if fx.touched[i]
      l[i] = return_cache(k,u.array[i],fx.array[i])
    end
  end
  ArrayBlock(g,fx.touched),l
end

function evaluate!(cache,k::LinearCombinationMap,u::ArrayBlock,fx::ArrayBlock)
  g,l = cache
  @check g.touched == fx.touched
  for i in eachindex(fx.array)
    if fx.touched[i]
      g.array[i] = evaluate!(l[i],k,u.array[i],fx.array[i])
    end
  end
  g
end

function Base.transpose(f::ArrayBlock{A,1} where A)
  fi = testitem(f)
  fit = transpose(fi)
  g = Matrix{typeof(fit)}(undef,(1,length(f.touched)))
  for i in eachindex(f.touched)
    if f.touched[i]
      g[i] = transpose(f.array[i])
    end
  end
  ArrayBlock(g,collect(transpose(f.touched)))
end

function Base.transpose(f::ArrayBlock{A,2} where A)
  fi = testitem(f)
  fit = transpose(fi)
  ni,nj = size(f)
  g = Matrix{typeof(fit)}(undef,(nj,ni))
  for i in 1:ni
    for j in 1:nj
      if f.touched[i,j]
        g[j,i] = transpose(f.array[i,j])
      end
    end
  end
  ArrayBlock(g,collect(transpose(f.touched)))
end

function return_cache(k::TransposeMap,f::ArrayBlock{A,1} where A)
  fi = testitem(f)
  li = return_cache(k,fi)
  fix = evaluate!(li,k,fi)
  l = Matrix{typeof(li)}(undef,(1,length(f.array)))
  g = Matrix{typeof(fix)}(undef,(1,length(f.array)))
  for i in eachindex(f.array)
    if f.touched[i]
      l[i] = return_cache(k,f.array[i])
    end
  end
  ArrayBlock(g,collect(transpose(f.touched))),l
end

function evaluate!(cache,k::TransposeMap,f::ArrayBlock{A,1} where A)
  g,l = cache
  @check g.touched == transpose(f.touched)
  for i in eachindex(f.array)
    if f.touched[i]
      g.array[i] = evaluate!(l[i],k,f.array[i])
    end
  end
  g
end

function integrate(f::ArrayBlock{A,N} where A,args...) where N
  fi = testitem(f)
  intfi = integrate(fi,args...)
  g = Array{typeof(intfi),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      g[i] = integrate(f.array[i],args...)
    end
  end
  ArrayBlock(g,f.touched)
end

function return_value(
  k::IntegrationMap,fx::ArrayBlock{A,N} where A,args...) where N
  fxi = testitem(fx)
  ufxi = return_value(k,fxi,args...)
  g = Array{typeof(ufxi),N}(undef,size(fx.array))
  for i in eachindex(fx.array)
    if fx.touched[i]
      g[i] = return_value(k,fx.array[i],args...)
    end
  end
  ArrayBlock(g,fx.touched)
end

function return_cache(
  k::IntegrationMap,fx::ArrayBlock{A,N} where A,args...) where N
  fxi = testitem(fx)
  li = return_cache(k,fxi,args...)
  ufxi = evaluate!(li,k,fxi,args...)
  l = Array{typeof(li),N}(undef,size(fx.array))
  g = Array{typeof(ufxi),N}(undef,size(fx.array))
  for i in eachindex(fx.array)
    if fx.touched[i]
      l[i] = return_cache(k,fx.array[i],args...)
    end
  end
  ArrayBlock(g,fx.touched),l
end

function evaluate!(cache,k::IntegrationMap,fx::ArrayBlock,args...)
  g,l = cache
  @check g.touched == fx.touched
  for i in eachindex(fx.array)
    if fx.touched[i]
      g.array[i] = evaluate!(l[i],k,fx.array[i],args...)
    end
  end
  g
end

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

function return_value(k::Broadcasting{typeof(∘)},f::ArrayBlock{A,N},h::Field) where {A,N}
  fi = testitem(f)
  fix = return_value(k,fi,h)
  g = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      g[i] = return_value(k,f.array[i],h)
    end
  end
  ArrayBlock(g,f.touched)
end

function return_cache(k::Broadcasting{typeof(∘)},f::ArrayBlock{A,N},h::Field) where {A,N}
  fi = testitem(f)
  li = return_cache(k,fi,h)
  fix = evaluate!(li,k,fi,h)
  l = Array{typeof(li),N}(undef,size(f.array))
  g = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      l[i] = return_cache(k,f.array[i],h)
    end
  end
  ArrayBlock(g,f.touched),l
end

function evaluate!(cache,k::Broadcasting{typeof(∘)},f::ArrayBlock,h::Field)
  g,l = cache
  @check g.touched == f.touched
  for i in eachindex(f.array)
    if f.touched[i]
      g.array[i] = evaluate!(l[i],k,f.array[i],h)
    end
  end
  g
end

function return_value(k::Broadcasting{<:Operation},f::ArrayBlock{A,N},h::Field) where {A,N}
  fi = testitem(f)
  fix = return_value(k,fi,h)
  g = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      g[i] = return_value(k,f.array[i],h)
    end
  end
  ArrayBlock(g,f.touched)
end

function return_cache(k::Broadcasting{<:Operation},f::ArrayBlock{A,N},h::Field) where {A,N}
  fi = testitem(f)
  li = return_cache(k,fi,h)
  fix = evaluate!(li,k,fi,h)
  l = Array{typeof(li),N}(undef,size(f.array))
  g = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      l[i] = return_cache(k,f.array[i],h)
    end
  end
  ArrayBlock(g,f.touched),l
end

function evaluate!(cache,k::Broadcasting{<:Operation},f::ArrayBlock,h::Field)
  g,l = cache
  @check g.touched == f.touched
  for i in eachindex(f.array)
    if f.touched[i]
      g.array[i] = evaluate!(l[i],k,f.array[i],h)
    end
  end
  g
end

function return_value(k::Broadcasting{<:Operation},h::Field,f::ArrayBlock{A,N}) where {A,N}
  fi = testitem(f)
  fix = return_value(k,h,fi)
  g = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      g[i] = return_value(k,h,f.array[i])
    end
  end
  ArrayBlock(g,f.touched)
end

function return_cache(k::Broadcasting{<:Operation},h::Field,f::ArrayBlock{A,N}) where {A,N}
  fi = testitem(f)
  li = return_cache(k,h,fi)
  fix = evaluate!(li,k,h,fi)
  l = Array{typeof(li),N}(undef,size(f.array))
  g = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      l[i] = return_cache(k,h,f.array[i])
    end
  end
  ArrayBlock(g,f.touched),l
end

function evaluate!(cache,k::Broadcasting{<:Operation},h::Field,f::ArrayBlock)
  g,l = cache
  @check g.touched == f.touched
  for i in eachindex(f.array)
    if f.touched[i]
      g.array[i] = evaluate!(l[i],k,h,f.array[i])
    end
  end
  g
end

function return_value(k::Broadcasting{<:Operation},h::ArrayBlock,f::ArrayBlock)
  evaluate(k,h,f)
end

function return_cache(k::Broadcasting{<:Operation},h::ArrayBlock,f::ArrayBlock{A,N}) where {A,N}
  @notimplemented
end

function evaluate!(cache,k::Broadcasting{<:Operation},h::ArrayBlock,f::ArrayBlock)
  @notimplemented
end

function return_value(k::BroadcastingFieldOpMap,f::ArrayBlock{A,N},g::AbstractArray) where {A,N}
  fi = testitem(f)
  fix = return_value(k,fi,g)
  h = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      h[i] = return_value(k,f.array[i],g)
    end
  end
  ArrayBlock(h,f.touched)
end

function return_cache(
  k::BroadcastingFieldOpMap,f::ArrayBlock{A,N},g::AbstractArray) where {A,N}
  fi = testitem(f)
  li = return_cache(k,fi,g)
  fix = evaluate!(li,k,fi,g)
  l = Array{typeof(li),N}(undef,size(f.array))
  h = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      l[i] = return_cache(k,f.array[i],g)
    end
  end
  ArrayBlock(h,f.touched),l
end

function evaluate!(cache,k::BroadcastingFieldOpMap,f::ArrayBlock,g::AbstractArray)
  h,l = cache
  @check h.touched == f.touched
  for i in eachindex(f.array)
    if f.touched[i]
      h.array[i] = evaluate!(l[i],k,f.array[i],g)
    end
  end
  h
end

function return_value(k::BroadcastingFieldOpMap,g::AbstractArray,f::ArrayBlock{A,N}) where {A,N}
  fi = testitem(f)
  fix = return_value(k,g,fi)
  h = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      h[i] = return_value(k,g,f.array[i])
    end
  end
  ArrayBlock(h,f.touched)
end

function return_cache(
  k::BroadcastingFieldOpMap,g::AbstractArray,f::ArrayBlock{A,N}) where {A,N}
  fi = testitem(f)
  li = return_cache(k,g,fi)
  fix = evaluate!(li,k,g,fi)
  l = Array{typeof(li),N}(undef,size(f.array))
  h = Array{typeof(fix),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      l[i] = return_cache(k,g,f.array[i])
    end
  end
  ArrayBlock(h,f.touched),l
end

function evaluate!(cache,k::BroadcastingFieldOpMap,g::AbstractArray,f::ArrayBlock)
  h,l = cache
  @check h.touched == f.touched
  for i in eachindex(f.array)
    if f.touched[i]
      h.array[i] = evaluate!(l[i],k,g,f.array[i])
    end
  end
  h
end

for op in (:+,:-,:*)
  @eval begin

    function return_value(k::Broadcasting{typeof($op)},f::ArrayBlock,g::ArrayBlock)
      return_value(BroadcastingFieldOpMap($op),f,g)
    end

    function return_cache(k::Broadcasting{typeof($op)},f::ArrayBlock,g::ArrayBlock)
      return_cache(BroadcastingFieldOpMap($op),f,g)
    end

    function evaluate!(cache,k::Broadcasting{typeof($op)},f::ArrayBlock,g::ArrayBlock)
      evaluate!(cache,BroadcastingFieldOpMap($op),f,g)
    end

  end
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
  evaluate(k,f,g)
end

function return_cache(k::Broadcasting{typeof(*)},f::ArrayBlock,g::Number)
  return_cache(k,g,f)
end

function evaluate!(cache,k::Broadcasting{typeof(*)},f::ArrayBlock,g::Number)
  evaluate!(cache,k,g,f)
end

function return_value(k::BroadcastingFieldOpMap,f::ArrayBlock,g::ArrayBlock)
  evaluate(k,f,g)
end

function return_cache(k::BroadcastingFieldOpMap,f::ArrayBlock,g::ArrayBlock)
  @notimplemented
end

function evaluate!(cache,k::BroadcastingFieldOpMap,f::ArrayBlock,g::ArrayBlock)
  @notimplemented
end

function return_value(
  k::BroadcastingFieldOpMap,f::ArrayBlock{A,N},g::ArrayBlock{B,N}) where {A,B,N}
  fi = testvalue(A)
  gi = testvalue(B)
  hi = return_value(k,fi,gi)
  a = Array{typeof(hi),N}(undef,size(f.array))
  fill!(a,hi)
  ArrayBlock(a,f.touched)
end

function return_cache(
  k::BroadcastingFieldOpMap,f::ArrayBlock{A,N},g::ArrayBlock{B,N}) where {A,B,N}
  @notimplementedif size(f) != size(g)
  fi = testvalue(A)
  gi = testvalue(B)
  ci = return_cache(k,fi,gi)
  hi = evaluate!(ci,k,fi,gi)
  m = ZeroBlockMap()
  a = Array{typeof(hi),N}(undef,size(f.array))
  b = Array{typeof(ci),N}(undef,size(f.array))
  zf = Array{typeof(return_cache(m,fi,gi))}(undef,size(f.array))
  zg = Array{typeof(return_cache(m,gi,fi))}(undef,size(f.array))
  t = map(|,f.touched,g.touched)
  for i in eachindex(f.array)
    if f.touched[i] && g.touched[i]
      b[i] = return_cache(k,f.array[i],g.array[i])
    elseif f.touched[i]
      _fi = f.array[i]
      zg[i] = return_cache(m,gi,_fi)
      _gi = evaluate!(zg[i],m,gi,_fi)
      b[i] = return_cache(k,_fi,_gi)
    elseif g.touched[i]
      _gi = g.array[i]
      zf[i] = return_cache(m,fi,_gi)
      _fi = evaluate!(zf[i],m,fi,_gi)
      b[i] = return_cache(k,_fi,_gi)
    end
  end
  ArrayBlock(a,t), b, zf, zg
end

function evaluate!(
  cache,k::BroadcastingFieldOpMap,f::ArrayBlock{A,N},g::ArrayBlock{B,N}) where {A,B,N}
  a,b, zf, zg = cache
  @check size(f) == size(g)
  @check size(a) == size(g)
  m = ZeroBlockMap()
  for i in eachindex(f.array)
    if f.touched[i] && g.touched[i]
      a.array[i] = evaluate!(b[i],k,f.array[i],g.array[i])
    elseif f.touched[i]
      fi = f.array[i]
      gi = evaluate!(zg[i],m,nothing,fi)
      a.array[i] = evaluate!(b[i],k,fi,gi)
    elseif g.touched[i]
      gi = g.array[i]
      fi = evaluate!(zf[i],m,nothing,gi)
      a.array[i] = evaluate!(b[i],k,fi,gi)
    end
  end
  a
end

function return_cache(
  k::BroadcastingFieldOpMap,f::ArrayBlock{A,1},g::ArrayBlock{B,2}) where {A,B}
  fi = testvalue(A)
  gi = testvalue(B)
  ci = return_cache(k,fi,gi)
  hi = evaluate!(ci,k,fi,gi)
  @check size(g.array,1) == 1 || size(g.array,2) == 0
  s = (size(f.array,1),size(g.array,2))
  a = Array{typeof(hi),2}(undef,s)
  b = Array{typeof(ci),2}(undef,s)
  t = fill(false,s)
  for j in 1:s[2]
    for i in 1:s[1]
      if f.touched[i] && g.touched[1,j]
        t[i,j] = true
        b[i,j] = return_cache(k,f.array[i],g.array[1,j])
      end
    end
  end
  ArrayBlock(a,t), b
end

function return_cache(
  k::BroadcastingFieldOpMap,f::ArrayBlock{A,2},g::ArrayBlock{B,1}) where {A,B}
  fi = testvalue(A)
  gi = testvalue(B)
  ci = return_cache(k,fi,gi)
  hi = evaluate!(ci,k,fi,gi)
  @check size(f.array,1) == 1 || size(f.array,2) == 0
  s = (size(g.array,1),size(f.array,2))
  a = Array{typeof(hi),2}(undef,s)
  b = Array{typeof(ci),2}(undef,s)
  t = fill(false,s)
  for j in 1:s[2]
    for i in 1:s[1]
      if f.touched[1,j] && g.touched[i]
        t[i,j] = true
        b[i,j] = return_cache(k,f.array[1,j],g.array[i])
      end
    end
  end
  ArrayBlock(a,t), b
end

function evaluate!(
  cache,k::BroadcastingFieldOpMap,f::ArrayBlock{A,1},g::ArrayBlock{B,2}) where {A,B}
  a,b = cache
  s = size(a.array)
  for j in 1:s[2]
    for i in 1:s[1]
      if f.touched[i] && g.touched[1,j]
        a.array[i,j] = evaluate!(b[i,j],k,f.array[i],g.array[1,j])
      end
    end
  end
  a
end

function evaluate!(
  cache,k::BroadcastingFieldOpMap,f::ArrayBlock{A,2},g::ArrayBlock{B,1}) where {A,B}
  a,b = cache
  s = size(a.array)
  for j in 1:s[2]
    for i in 1:s[1]
      if f.touched[1,j] && g.touched[i]
        a.array[i,j] = evaluate!(b[i,j],k,f.array[1,j],g.array[i])
      end
    end
  end
  a
end

function return_value(
  k::BroadcastingFieldOpMap,a::(ArrayBlock{A,N} where A)...) where N
  evaluate(k,a...)
end

function return_cache(
  k::BroadcastingFieldOpMap,a::(ArrayBlock{A,N} where A)...) where N
  a1 = first(a)
  @notimplementedif any(ai->size(ai)!=size(a1),a)
  ais = map(ai->testvalue(eltype(ai)),a)
  ci = return_cache(k,ais...)
  bi = evaluate!(ci,k,ais...)
  c = Array{typeof(ci),N}(undef,size(a1))
  array = Array{typeof(bi),N}(undef,size(a1))
  for i in eachindex(a1.array)
    @notimplementedif any(ai->ai.touched[i]!=a1.touched[i],a)
    if a1.touched[i]
      _ais = map(ai->ai.array[i],a)
      c[i] = return_cache(k,_ais...)
    end
  end
  ArrayBlock(array,a1.touched), c
end

function evaluate!(
  cache,k::BroadcastingFieldOpMap,a::(ArrayBlock{A,N} where A)...) where N
  a1 = first(a)
  @notimplementedif any(ai->size(ai)!=size(a1),a)
  r,c = cache
  for i in eachindex(a1.array)
    @notimplementedif any(ai->ai.touched[i]!=a1.touched[i],a)
    if a1.touched[i]
      ais = map(ai->ai.array[i],a)
      r.array[i] = evaluate!(c[i],k,ais...)
    end
  end
  r
end

function return_value(k::BroadcastingFieldOpMap,a::ArrayBlock...)
  evaluate(k,a...)
end

function return_cache(k::BroadcastingFieldOpMap,a::ArrayBlock...)
  @notimplemented
end

function evaluate!(cache,k::BroadcastingFieldOpMap,a::ArrayBlock...)
  @notimplemented
end

function return_value(
  k::BroadcastingFieldOpMap,a::Union{ArrayBlock,AbstractArray}...)
  evaluate(k,a...)
end

function return_cache(
  k::BroadcastingFieldOpMap,a::Union{ArrayBlock,AbstractArray}...)

  function _replace_nz_blocks(a::ArrayBlock,bi::AbstractArray)
    N = ndims(a.array)
    array = Array{typeof(bi),N}(undef,size(a))
    for i in eachindex(a.array)
      if a.touched[i]
        array[i] = bi
      end
    end
    ArrayBlock(array,a.touched)
  end

  function _replace_nz_blocks(a::ArrayBlock,bi::ArrayBlock)
    bi
  end

  inds = findall(ai->isa(ai,ArrayBlock),a)
  @notimplementedif length(inds) == 0
  a1 = a[inds[1]]
  b = map(ai->_replace_nz_blocks(a1,ai),a)
  c = return_cache(k,b...)
  c,b
end

function evaluate!(
  cache,k::BroadcastingFieldOpMap,a::Union{ArrayBlock,AbstractArray}...)

  function _replace_nz_blocks!(a::ArrayBlock,bi::AbstractArray)
    for i in eachindex(a.array)
      if a.touched[i]
        a.array[i] = bi
      end
    end
    a
  end

  function _replace_nz_blocks!(a::ArrayBlock,bi::ArrayBlock)
    bi
  end

  c, b = cache
  d = map((i,j)->_replace_nz_blocks!(i,j),b,a)
  evaluate!(c,k,d...)
end

function Base.:+(a::ArrayBlock,b::ArrayBlock)
  BroadcastingFieldOpMap(+)(a,b)
end

function Base.:-(a::ArrayBlock,b::ArrayBlock)
  BroadcastingFieldOpMap(-)(a,b)
end

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

function Base.:*(a::ArrayBlock,b::Number)
  b*a
end

function return_value(::typeof(*),a::ArrayBlock{A,2},b::ArrayBlock{B,1}) where {A,B}
  ai = testvalue(A)
  bi = testvalue(B)
  ri = return_value(*,ai,bi)
  array = Vector{typeof(ri)}(undef,size(a.array,1))
  touched = fill(false,size(a.array,1))
  ArrayBlock(array,touched)
end

function Base.:*(a::ArrayBlock{A,2},b::ArrayBlock{B,1}) where {A,B}
  @check size(a.array,2) == size(b.array,1)
  ai = testvalue(A)
  bi = testvalue(B)
  ri = ai*bi
  array = Vector{typeof(ri)}(undef,size(a.array,1))
  touched = fill(false,size(a.array,1))
  ni,nj = size(a.array)
  for i in 1:ni
    for j in 1:nj
      if a.touched[i,j] && b.touched[j]
        if !touched[i]
          array[i] = a.array[i,j]*b.array[j]
          touched[i] = true
        else
          array[i] = array[i] + a.array[i,j]*b.array[j]
        end
      end
    end
  end
  ArrayBlock(array,touched)
end

# Mostly used for v*tranpose(v)
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
  ai = testvalue(A)
  bi = testvalue(B)
  ri = ai*bi
  ni = size(a.array,1)
  nj = size(b.array,2)
  nk = size(a.array,2)
  s = (ni,nj)
  array = Matrix{typeof(ri)}(undef,s)
  touched = fill(false,s)
  ni,nj = size(a.array)
  for i in 1:ni
    for j in 1:nj
      for k in 1:nk
        if a.touched[i,k] && b.touched[k,j]
          if !touched[i,j]
            array[i,j] = a.array[i,k]*b.array[k,j]
            touched[i,j] = true
          else
            array[i,j] = array[i,j] + a.array[i,k]*b.array[k,j]
          end
        end
      end
    end
  end
  ArrayBlock(array,touched)
end

function LinearAlgebra.rmul!(a::ArrayBlock,β)
  for i in eachindex(a.touched)
    if a.touched[i]
      rmul!(a.array[i],β)
    end
  end
end

function _zero_entries!(a::AbstractArray)
  fill!(a,zero(eltype(a)))
end

function _zero_entries!(a::ArrayBlock)
  for i in eachindex(a.touched)
    if a.touched[i]
      _zero_entries!(a.array[i])
    end
  end
end

function LinearAlgebra.mul!(c::ArrayBlock,a::ArrayBlock,b::ArrayBlock)
  _zero_entries!(c)
  mul!(c,a,b,1,0)
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
          _mymul!(cij,aik,bkj,α)
          # hack to avoid allocations
          #mul!(cij,aik,bkj,α,1)
        end
      end
    end
  end
  c
end

function  _setsize_mul!(c,a::AbstractMatrix,b::AbstractVector)
  ni = size(a,1)
  @check length(b) == size(a,2)
  setsize!(c,(ni,))
  c
end

function  _setsize_mul!(c,a::AbstractMatrix,b::AbstractMatrix)
  ni = size(a,1)
  nj = size(b,2)
  @check size(a,2) == size(b,1)
  setsize!(c,(ni,nj))
  c
end

function  _setsize_mul!(c,a::MatrixBlock,b::VectorBlock)
  ni,nj = size(a)
  @check length(b) == nj
  @check length(c) == ni
  for i in 1:ni
    for j in 1:nj
      if a.touched[i,j] && b.touched[j]
        _setsize_mul!(c.array[i],a.array[i,j],b.array[j])
      end
    end
  end
  c
end

function  _setsize_mul!(c,a::MatrixBlock,b::MatrixBlock)
  ni,nk = size(a)
  nk2,nj = size(b)
  @check nk == nk2
  @check size(c) == (ni,nj)
  for i in 1:ni
    for j in 1:nj
      for k in 1:nk
        if a.touched[i,k] && b.touched[k,j]
          _setsize_mul!(c.array[i,j],a.array[i,k],b.array[k,j])
        end
      end
    end
  end
  c
end

function return_cache(::typeof(*),a::ArrayBlock,b::ArrayBlock)
  c1 = CachedArray(a*b)
  c2 = return_cache(unwrap_cached_array,c1)
  (c1,c2)
end

function evaluate!(cache,::typeof(*),a::ArrayBlock,b::ArrayBlock)
  c1,c2 = cache
  _setsize_mul!(c1,a,b)
  c = evaluate!(c2,unwrap_cached_array,c1)
  mul!(c,a,b)
  c
end

function  _setsize_as!(d,a::AbstractArray)
  setsize!(d,size(a))
end

function  _setsize_as!(d,a::ArrayBlock)
  @check size(d) == size(a)
  for i in eachindex(a.array)
    if a.touched[i]
      _setsize_as!(d.array[i],a.array[i])
    end
  end
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
  _setsize_as!(c1,c)
  _setsize_mul!(c1,a,b)
  d = evaluate!(c2,unwrap_cached_array,c1)
  copyto!(d,c)
  iszero(k.α) && isone(k.β) && return d
  mul!(d,a,b,k.α,k.β)
  d
end

function Base.copyto!(d::ArrayBlock,c::ArrayBlock)
  _zero_entries!(d)
  for i in eachindex(c.array)
    if c.touched[i]
      copyto!(d.array[i],c.array[i])
    end
  end
  d
end

# hack to avoid allocations
function _mymul!(
  c::ArrayBlock{C,2} where C,
  a::ArrayBlock{A,2} where A,
  b::ArrayBlock{B,2} where B,
  α::Number)
  mul!(c,a,b,α,1)
end

function _mymul!(cIJ,aIK,bKJ,α)
  @boundscheck begin
    @assert size(cIJ,1) == size(aIK,1)
    @assert size(cIJ,2) == size(bKJ,2)
    @assert size(aIK,2) == size(bKJ,1)
  end
  for i in 1:size(cIJ,1)
    for j in 1:size(bKJ,2)
      for k in 1:size(bKJ,1)
        @inbounds cIJ[i,j] += α*aIK[i,k]*bKJ[k,j]
      end
    end
  end
end

# Assembly related

for T in (:AddEntriesMap,:TouchEntriesMap)
  @eval begin

    function return_cache(
      k::$T,A,v::MatrixBlock,I::VectorBlock,J::VectorBlock)

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

    function evaluate!(
      cache, k::$T,A,v::MatrixBlock,I::VectorBlock,J::VectorBlock)
      ni,nj = size(v.touched)
      for j in 1:nj
        for i in 1:ni
          if v.touched[i,j]
            evaluate!(cache[i,j],k,A,v.array[i,j],I.array[i],J.array[j])
          end
        end
      end
    end

    function return_cache(
      k::$T,A,v::VectorBlock,I::VectorBlock)

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

    function evaluate!(
      cache, k::$T,A,v::VectorBlock,I::VectorBlock)
      ni = length(v.touched)
      for i in 1:ni
        if v.touched[i]
          evaluate!(cache[i],k,A,v.array[i],I.array[i])
        end
      end
    end

  end
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

_zero_entries!(a::ArrayBlockView) = _zero_entries!(a.array)

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

# Autodiff - Skeleton + SingleField

function return_cache(k::Arrays.AutoDiffMap,cfg::ForwardDiff.JacobianConfig,ydual::VectorBlock)
  i = findfirst(ydual.touched)
  yi = ydual.array[i]
  ci = return_cache(k,cfg,yi)
  ri = evaluate!(ci,k,cfg,yi)
  cache = Vector{typeof(ci)}(undef,length(ydual.array))
  array = Vector{typeof(ri)}(undef,length(ydual.array))
  for i in eachindex(ydual.array)
    if ydual.touched[i]
      cache[i] = return_cache(k,cfg,ydual.array[i])
    end
  end
  result = ArrayBlock(array,ydual.touched)
  return result, cache
end

function evaluate!(cache,k::Arrays.AutoDiffMap,cfg::ForwardDiff.JacobianConfig,ydual::VectorBlock)
  r, c = cache
  for i in eachindex(ydual.array)
    if ydual.touched[i]
      r.array[i] = evaluate!(c[i],k,cfg,ydual.array[i])
    end
  end
  return r
end

# Autodiff - MultiField

struct BlockConfig{C,T,V,N,D,O} <: ForwardDiff.AbstractConfig{N}
  seeds::NTuple{N,ForwardDiff.Partials{N,V}}
  duals::D
  offsets::O
end

function BlockConfig(
  op::Union{typeof(ForwardDiff.gradient),typeof(ForwardDiff.jacobian)},
  f::F,
  x::Union{VectorBlock{<:AbstractArray{V}},VectorBlock{<:VectorBlock{<:AbstractArray{V}}}},
  ::T = ForwardDiff.Tag(f,V)
) where {F,V,T}
  offsets, N = block_offsets(x, 0)
  seeds = ForwardDiff.construct_seeds(ForwardDiff.Partials{N,V})
  duals = similar(x, ForwardDiff.Dual{T,V,N})
  BlockConfig{typeof(op),T,V,N,typeof(duals),typeof(offsets)}(seeds,duals,offsets)
end

@inline block_offsets(x::Vector, offset) = offset, offset + length(x)

function block_offsets(x::VectorBlock{A}, offset) where A
  offsets = ()
  for i in eachindex(x.touched)
    if x.touched[i]
      @inbounds offsets_i, offset = block_offsets(x.array[i], offset)
    else
      offsets_i = -1
    end
    offsets = (offsets...,offsets_i)
  end
  return offsets, offset
end

for F in (ForwardDiff.gradient,ForwardDiff.jacobian)
  @eval begin
    function Arrays.return_cache(k::ConfigMap{typeof($F)},x::VectorBlock)
      return BlockConfig($F,k.tag,x)
    end
  end
end

function evaluate!(cache,k::DualizeMap,cfg::BlockConfig,x)
  xdual, seeds, offsets = cfg.duals, cfg.seeds, cfg.offsets
  seed_block!(xdual, x, seeds, offsets)
  return xdual
end

function Arrays.return_cache(::AutoDiffMap,cfg::BlockConfig{typeof(ForwardDiff.gradient),T},ydual) where T
  ydual isa Real || throw(ForwardDiff.GRAD_ERROR)
  result = similar(cfg.duals, ForwardDiff.valtype(ydual))
  return result
end

function Arrays.evaluate!(result,::AutoDiffMap,cfg::BlockConfig{typeof(ForwardDiff.gradient),T},ydual) where T
  extract_gradient_block!(T, result, ydual, cfg.offsets)
  return result
end

function return_cache(::AutoDiffMap,cfg::BlockConfig{typeof(ForwardDiff.jacobian),T},ydual) where T
  ydual isa VectorBlock || throw(ForwardDiff.JACOBIAN_ERROR)
  result = _alloc_jacobian(ydual,cfg.duals)
  return result
end

function evaluate!(result,::AutoDiffMap,cfg::BlockConfig{typeof(ForwardDiff.jacobian),T},ydual) where T
  extract_jacobian_block!(T, result, ydual, cfg.offsets)
  return result
end

function _alloc_jacobian(ydual::Vector,xdual::Vector)
  T = ForwardDiff.valtype(eltype(ydual))
  zeros(T,length(ydual),length(xdual))
end

# Skeleton + Multifield: The VectorBlock corresponds to +/-
function _alloc_jacobian(ydual::VectorBlock,xdual::Vector)
  i = findfirst(ydual.touched)
  ai = _alloc_jacobian(ydual.array[i],xdual)
  ni = size(ydual.array,1)
  array = Vector{typeof(ai)}(undef,ni)
  for i in 1:ni
    if ydual.touched[i]
      array[i] = _alloc_jacobian(ydual.array[i],xdual)
    end
  end
  ArrayBlock(array,ydual.touched)
end

function _alloc_jacobian(ydual::VectorBlock,xdual::VectorBlock)
  i = findfirst(ydual.touched)
  j = findfirst(xdual.touched)
  ai = _alloc_jacobian(ydual.array[i],xdual.array[j])

  ni, nj = size(ydual.array,1), size(xdual.array,1)
  array = Matrix{typeof(ai)}(undef,ni,nj)
  touched = fill(false,ni,nj)
  for i in 1:ni
    for j in 1:nj
      if ydual.touched[i] && xdual.touched[j]
        array[i,j]   = _alloc_jacobian(ydual.array[i],xdual.array[j])
        touched[i,j] = true
      end
    end
  end
  ArrayBlock(array,touched)
end

function seed_block!(
  duals::VectorBlock{A}, x::VectorBlock{B}, seeds::NTuple{N,ForwardDiff.Partials{N}}, offsets
) where {N,A,B}
  for i in eachindex(duals.touched)
    if duals.touched[i]
      @check x.touched[i]
      @inbounds seed_block!(duals.array[i], x.array[i], seeds, offsets[i])
    end
  end
  return duals
end

function seed_block!(
  duals::AbstractArray{ForwardDiff.Dual{T,V,N}}, x, seeds::NTuple{N,ForwardDiff.Partials{N,V}}, offset
) where {T,V,N}
  for j in eachindex(duals)
    @inbounds duals[j] = ForwardDiff.Dual{T,V,N}(x[j], seeds[j+offset])
  end
  return duals
end

function extract_gradient_block!(::Type{T}, result::VectorBlock{A}, dual, offsets) where {T,A}
  for i in eachindex(result.touched)
    if result.touched[i]
      @inbounds extract_gradient_block!(T, result.array[i], dual, offsets[i])
    end
  end
  return result
end

function extract_gradient_block!(::Type{T}, result::AbstractArray, dual::ForwardDiff.Dual, offset) where {T}
  for j in eachindex(result)
    @inbounds result[j] = ForwardDiff.partials(T,dual,j+offset)
  end
  return result
end

function extract_gradient_block!(::Type{T}, result::AbstractArray, dual::Real, offset) where {T}
  fill!(result,zero(dual))
  return result
end

function extract_jacobian_block!(::Type{T}, result::MatrixBlock{A}, dual::VectorBlock{B}, offsets) where {T,A,B}
  for i in axes(result.touched,1)
    for j in axes(result.touched,2)
      if result.touched[i,j]
        @inbounds extract_jacobian_block!(T, result.array[i,j], dual.array[i], offsets[j])
      end
    end
  end
  return result
end

# Skeleton + Multifield: The VectorBlocks correspond to +/-
function extract_jacobian_block!(::Type{T}, result::VectorBlock{A}, dual::VectorBlock{B}, offset) where {T,A,B}
  for i in axes(result.touched,1)
    if result.touched[i]
      @check dual.touched[i]
      @inbounds extract_jacobian_block!(T, result.array[i], dual.array[i], offset)
    end
  end
  return result
end

function extract_jacobian_block!(::Type{T}, result::AbstractArray, dual::AbstractArray{<:ForwardDiff.Dual}, offset) where {T}
  for k in axes(result,1)
    for l in axes(result,2)
      @inbounds result[k,l] = ForwardDiff.partials(T,dual[k],l+offset)
    end
  end
  return result
end

function extract_jacobian_block!(::Type{T}, result::AbstractArray, dual::AbstractArray{<:Real}, offset) where {T}
  fill!(result,zero(eltype(dual)))
  return result
end
