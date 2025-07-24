
# Optimisations for lazy_map with BlockMap

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

# evaluate + ArrayBlock

function return_value(f::ArrayBlock{A,N},x) where {A,N}
  fi = testitem(f)
  fix = return_value(fi,x)
  g = Array{typeof(fix),N}(undef,size(f.array))
  fill!(g,fix)
  ArrayBlock(g,f.touched)
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

# linear_combination + ArrayBlock

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

function return_value(k::LinearCombinationMap,u::ArrayBlock,fx::ArrayBlock)
  i::Int = findfirst(fx.touched)
  fxi = fx.array[i]
  ui = u.array[i]
  ufxi = return_value(k,ui,fxi)
  g = Vector{typeof(ufxi)}(undef,size(fx.array))
  fill!(g,ufxi)
  ArrayBlock(g,fx.touched)
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

# transpose + ArrayBlock

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

function return_value(k::TransposeMap,f::ArrayBlock{A,1} where A)
  fi = testitem(f)
  fix = return_value(k,fi)
  g = Matrix{typeof(fix)}(undef,(1,length(f.array)))
  fill!(g,fix)
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

# integrate + ArrayBlock

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

# Broadcasting + ArrayBlock for composition of fields

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

# Show error for everything else
#
# This is not implemented for a good reason, as it is not clear how to return a 
# concrete type that would make sense.
#
# For example, let's assume we have a MultiFieldFEBasis with components xh = (uh, ph)
# and we want to have yh = uh + ph
# As long as both evaluate to the same type of scalar T, we can evaluate uh(x) and ph(x)
# separately and then add them to obtain yh(x). The resulting ArrayBlock would have 
# a concrete block type Array{T,N}. 
# In general, the final evaluated weakform has to be of the same scalar type regardless 
# of the basis functions used in each variable, since it needs to be integrated. 
# So the final result can always be broadcasted to a concrete block type that we can wrap in 
# an ArrayBlock.
# However, this is NOT true for the fields themselves. If I do yh = uh + ph, what I should
# return is yh = [uh, 0] + [0, ph] = [uh, ph]. Therefore in each block I would have a
# different field type, which would then make the ArrayBlock type ambiguous.
#
# This is why all of these cases are disabled, and one should be careful when defining new stuff.

function return_value(k::Broadcasting{<:Operation},h::ArrayBlock,f::ArrayBlock)
  evaluate(k,h,f)
end

function return_cache(k::Broadcasting{<:Operation},h::ArrayBlock,f::ArrayBlock{A,N}) where {A,N}
  @notimplemented
end

function evaluate!(cache,k::Broadcasting{<:Operation},h::ArrayBlock,f::ArrayBlock)
  @notimplemented 
end

# BroadcastingFieldOpMap for ArrayBlocks

## Case 1: ArrayBlock + AbstractArray

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

function return_cache(k::BroadcastingFieldOpMap,f::ArrayBlock{A,N},g::AbstractArray) where {A,N}
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

function return_cache(k::BroadcastingFieldOpMap,g::AbstractArray,f::ArrayBlock{A,N}) where {A,N}
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

# Case 2: ArrayBlock + ArrayBlock of same dimension

function return_value(k::BroadcastingFieldOpMap,f::ArrayBlock{A,N},g::ArrayBlock{B,N}) where {A,B,N}
  fi = testvalue(A)
  gi = testvalue(B)
  hi = return_value(k,fi,gi)
  a = Array{typeof(hi),N}(undef,size(f.array))
  fill!(a,hi)
  ArrayBlock(a,f.touched)
end

function return_cache(k::BroadcastingFieldOpMap,f::ArrayBlock{A,N},g::ArrayBlock{B,N}) where {A,B,N}
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

function evaluate!(cache,k::BroadcastingFieldOpMap,f::ArrayBlock{A,N},g::ArrayBlock{B,N}) where {A,B,N}
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

# Case 3: ArrayBlock + ArrayBlock of different dimensions
# There are 2 subcases: 
#   - ArrayBlock{A,1} vs ArrayBlock{B,2} comes from test vs trial
#   - ArrayBlock{A,2} vs ArrayBlock{B,1} comes from trial vs test

function return_cache(k::BroadcastingFieldOpMap,f::ArrayBlock{A,1},g::ArrayBlock{B,2}) where {A,B}
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

function return_cache(k::BroadcastingFieldOpMap,f::ArrayBlock{A,2},g::ArrayBlock{B,1}) where {A,B}
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

function evaluate!(cache,k::BroadcastingFieldOpMap,f::ArrayBlock{A,1},g::ArrayBlock{B,2}) where {A,B}
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

function evaluate!(cache,k::BroadcastingFieldOpMap,f::ArrayBlock{A,2},g::ArrayBlock{B,1}) where {A,B}
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

# Case 4: Catch all operations between ArrayBlocks of same dimension

function return_value(k::BroadcastingFieldOpMap,a::(ArrayBlock{A,N} where A)...) where N
  evaluate(k,a...)
end

function return_cache(k::BroadcastingFieldOpMap,a::(ArrayBlock{A,N} where A)...) where N
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

function evaluate!(cache,k::BroadcastingFieldOpMap,a::(ArrayBlock{A,N} where A)...) where N
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

# Case 5: Catch-all version, where any combination of ArrayBlock and AbstractArray is allowed.
# However: We have noticed some allocations in this case. 

function return_value(k::BroadcastingFieldOpMap,a::Union{ArrayBlock,AbstractArray}...)
  evaluate(k,a...)
end

function return_cache(k::BroadcastingFieldOpMap,a::Union{ArrayBlock,AbstractArray}...)
  i = findfirst(Base.Fix2(isa,ArrayBlock),a)
  @notimplementedif isnothing(i)
  m = Arrays.MatchingBlockMap(a[i])
  cm = map(ai -> return_cache(m, ai), a)
  b = map((ci,ai) -> evaluate!(ci, m, ai), cm, a)
  c = return_cache(k,b...)
  return c, m, cm
end

function evaluate!(cache,k::BroadcastingFieldOpMap,a::Union{ArrayBlock,AbstractArray}...)
  c, m, cm = cache
  b = map((ci,ai) -> evaluate!(ci, m, ai), cm, a)
  evaluate!(c,k,b...)
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

function Base.:+(a::ArrayBlock,b::ArrayBlock)
  BroadcastingFieldOpMap(+)(a,b)
end

function Base.:-(a::ArrayBlock,b::ArrayBlock)
  BroadcastingFieldOpMap(-)(a,b)
end

function Base.:-(a::ArrayBlock)
  BroadcastingFieldOpMap(-)(a)
end
