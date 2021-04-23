
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


# TODO rename it to ArrayBlock when we don't need to depend on BlockArrays
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

Base.size(b::ArrayBlock) = size(b.array)
Base.length(b::ArrayBlock) = length(b.array)
Base.eltype(::Type{<:ArrayBlock{A}}) where A = A
Base.eltype(b::ArrayBlock{A}) where A = A
Base.ndims(b::ArrayBlock{A,N}) where {A,N} = N
Base.ndims(::Type{ArrayBlock{A,N}}) where {A,N} = N
@inline function Base.getindex(b::ArrayBlock,i...)
  if !b.touched[i...]
    return nothing
  end
  b.array[i...]
end
#@inline function Base.setindex!(b::ArrayBlock,v,i...)
#  @check b.touched[i...] "Only touched entries can be set"
#  b.array[i...] = v
#end
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

#LinearAlgebra.promote_leaf_eltypes(a::ArrayBlock) = LinearAlgebra.promote_leaf_eltypes(a.array)

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

struct BlockMap{N}
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

function  return_cache(k::TransposeMap,f::ArrayBlock{A,1} where A)
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

function  evaluate!(cache,k::TransposeMap,f::ArrayBlock)
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

function return_value(k::Broadcasting,f::ArrayBlock)
  evaluate(k,f)
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

function return_value(k::Broadcasting{typeof(∘)},f::ArrayBlock,h::Field)
  evaluate(k,f,h)
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

function return_value(k::Broadcasting{<:Operation},f::ArrayBlock,h::Field)
  evaluate(k,f,h)
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

function return_value(k::Broadcasting{<:Operation},h::Field,f::ArrayBlock)
  evaluate(k,h,f)
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

function return_value(k::BroadcastingFieldOpMap,f::ArrayBlock,g::AbstractArray)
  evaluate(k,f,g)
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

function return_value(k::BroadcastingFieldOpMap,g::AbstractArray,f::ArrayBlock)
  evaluate(k,g,f)
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
      evaluate(k,f,g)
    end
    
    function return_cache(k::Broadcasting{typeof($op)},f::ArrayBlock,g::ArrayBlock)
      return_cache(BroadcastingFieldOpMap($op),f,g)
    end
    
    function evaluate!(cache,k::Broadcasting{typeof($op)},f::ArrayBlock,g::ArrayBlock)
      evaluate!(cache,BroadcastingFieldOpMap($op),f,g)
    end

  end
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

function Base.:+(a::ArrayBlock,b::ArrayBlock)
  BroadcastingFieldOpMap(+)(a,b)
end

function Base.:-(a::ArrayBlock,b::ArrayBlock)
  BroadcastingFieldOpMap(-)(a,b)
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

function Algebra.scale_entries!(a::ArrayBlock,β)
  for i in eachindex(a.touched)
    if a.touched[i]
      scale_entries!(a.array[i],β)
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
      scale_entries!(c.array[i],β)
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
        scale_entries!(c.array[i,j],β)
      end
      for k in 1:nk
        if a.touched[i,k] && b.touched[k,j]
          cij = c.array[i,j]
          aik = a.array[i,k]
          bkj = b.array[k,j]
          Arrays.mymul!(cij,aik,bkj,α)
          # hack to avoid allocations
          #mul!(cij,aik,bkj,α,1)
        end
      end
    end
  end
  c
end

# hack to avoid allocations
function Arrays.mymul!(
  c::ArrayBlock{C,2} where C,
  a::ArrayBlock{A,2} where A,
  b::ArrayBlock{B,2} where B,
  α::Number)
  mul!(c,a,b,α,1)
end

function return_cache(::typeof(*),a::ArrayBlock,b::ArrayBlock)
  a*b
end

function evaluate!(c,::typeof(*),a::ArrayBlock,b::ArrayBlock)
  mul!(c,a,b)
  c
end

function return_cache(k::MulAddMap,a::ArrayBlock,b::ArrayBlock,c::ArrayBlock)
  d = a*b+c
  d
end

function evaluate!(d,k::MulAddMap,a::ArrayBlock,b::ArrayBlock,c::ArrayBlock)
  copyto!(d,c)
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


