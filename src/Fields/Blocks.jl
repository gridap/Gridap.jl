
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

# TODO rename it to GBlock when we don't need to depend on BlockArrays
struct GBlock{A,N}
  array::Array{A,N}
  touched::Array{Bool,N}
  function GBlock(array::Array{A,N},touched::Array{Bool,N}) where {A,N}
    @check size(array) == size(touched)
    new{A,N}(array,touched)
  end
end

Base.size(b::GBlock) = size(b.array)
Base.length(b::GBlock) = length(b.array)
Base.eltype(b::GBlock{A}) where A = A
Base.ndims(b::GBlock{A,N}) where {A,N} = N
@inline function Base.getindex(b::GBlock,i...)
  if !b.touched[i...]
    return nothing
  end
  b.array[i...]
end
#@inline function Base.setindex!(b::GBlock,v,i...)
#  @check b.touched[i...] "Only touched entries can be set"
#  b.array[i...] = v
#end
function Base.show(io::IO,o::GBlock)
  print(io,"GBlock($(o.array), $(o.touched))")
end
#function Base.show(io::IO,k::MIME"text/plain",o::GBlock)
#  show(io,k, map((a,b)->b ? a : nothing,o.array,o.touched))
#end

function Arrays.testitem(f::GBlock{A,1}) where A
  @notimplementedif !isconcretetype(A)
  @notimplementedif length(f.array) == 0
  i::Int = findfirst(f.touched)
  fi = f.array[i]
  fi
end

function Arrays.testitem(f::GBlock{A,N}) where {A,N}
  @notimplementedif !isconcretetype(A)
  @notimplementedif length(f.array) == 0
  i::CartesianIndex{N} = findfirst(f.touched)
  fi = f.array[i]
  fi
end

function Arrays.testvalue(::Type{GBlock{A,N}}) where {A,N}
  s = ntuple(i->0,Val(N))
  array = Array{A,N}(undef,s)
  touched = Array{Bool,N}(undef,s)
  GBlock(array,touched)
end

function return_cache(f::GBlock{A,N},x) where {A,N}
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
  GBlock(g,f.touched),l
end

function evaluate!(cache,f::GBlock,x)
  g,l = cache
  @check g.touched == f.touched
  for i in eachindex(f.array)
    if f.touched[i]
      g.array[i] = evaluate!(l[i],f.array[i],x)
    end
  end
  g
end

function linear_combination(u::GBlock,f::GBlock)
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
  GBlock(g,f.touched)
end

function return_cache(k::LinearCombinationMap,u::GBlock,fx::GBlock)
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
  GBlock(g,fx.touched),l
end

function evaluate!(cache,k::LinearCombinationMap,u::GBlock,fx::GBlock)
  g,l = cache
  @check g.touched == fx.touched
  for i in eachindex(fx.array)
    if fx.touched[i]
      g.array[i] = evaluate!(l[i],k,u.array[i],fx.array[i])
    end
  end
  g
end

function Base.transpose(f::GBlock{A,1} where A)
  fi = testitem(f)
  fit = transpose(fi)
  g = Matrix{typeof(fit)}(undef,(1,length(f.touched)))
  for i in eachindex(f.touched)
    if f.touched[i]
      g[i] = transpose(f.array[i])
    end
  end
  GBlock(g,collect(transpose(f.touched)))
end

function  return_cache(k::TransposeMap,f::GBlock{A,1} where A)
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
  GBlock(g,collect(transpose(f.touched))),l
end

function  evaluate!(cache,k::TransposeMap,f::GBlock)
  g,l = cache
  @check g.touched == transpose(f.touched)
  for i in eachindex(f.array)
    if f.touched[i]
      g.array[i] = evaluate!(l[i],k,f.array[i])
    end
  end
  g
end

function integrate(f::GBlock{A,N} where A,args...) where N
  fi = testitem(f)
  intfi = integrate(fi,args...)
  g = Array{typeof(intfi),N}(undef,size(f.array))
  for i in eachindex(f.array)
    if f.touched[i]
      g[i] = integrate(f.array[i],args...)
    end
  end
  GBlock(g,f.touched)
end

function return_cache(
  k::IntegrationMap,fx::GBlock{A,N} where A,args...) where N
  i::Int = findfirst(fx.touched)
  fxi = fx.array[i]
  li = return_cache(k,fxi,args...)
  ufxi = evaluate!(li,k,fxi,args...)
  l = Array{typeof(li),N}(undef,size(fx.array))
  g = Array{typeof(ufxi),N}(undef,size(fx.array))
  for i in eachindex(fx.array)
    if fx.touched[i]
      l[i] = return_cache(k,fx.array[i],args...)
    end
  end
  GBlock(g,fx.touched),l
end

function evaluate!(cache,k::IntegrationMap,fx::GBlock,args...)
  g,l = cache
  @check g.touched == fx.touched
  for i in eachindex(fx.array)
    if fx.touched[i]
      g.array[i] = evaluate!(l[i],k,fx.array[i],args...)
    end
  end
  g
end

function return_value(k::Broadcasting,f::GBlock)
  evaluate(k,f)
end

function return_cache(k::Broadcasting,f::GBlock{A,N}) where {A,N}
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
  GBlock(g,f.touched),l
end

function evaluate!(cache,k::Broadcasting,f::GBlock)
  g,l = cache
  @check g.touched == f.touched
  for i in eachindex(f.array)
    if f.touched[i]
      g.array[i] = evaluate!(l[i],k,f.array[i])
    end
  end
  g
end

function return_value(k::Broadcasting{typeof(∘)},f::GBlock,h::Field)
  evaluate(k,f,h)
end

function return_cache(k::Broadcasting{typeof(∘)},f::GBlock{A,N},h::Field) where {A,N}
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
  GBlock(g,f.touched),l
end

function evaluate!(cache,k::Broadcasting{typeof(∘)},f::GBlock,h::Field)
  g,l = cache
  @check g.touched == f.touched
  for i in eachindex(f.array)
    if f.touched[i]
      g.array[i] = evaluate!(l[i],k,f.array[i],h)
    end
  end
  g
end

function Fields.return_value(k::BroadcastingFieldOpMap,f::GBlock,g::AbstractArray)
  evaluate(k,f,g)
end

function Fields.return_cache(
  k::BroadcastingFieldOpMap,f::GBlock{A,N},g::AbstractArray) where {A,N}
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
  GBlock(h,f.touched),l
end

function Fields.evaluate!(cache,k::BroadcastingFieldOpMap,f::GBlock,g::AbstractArray)
  h,l = cache
  @check h.touched == f.touched
  for i in eachindex(f.array)
    if f.touched[i]
      h.array[i] = evaluate!(l[i],k,f.array[i],g)
    end
  end
  h
end

function Fields.return_value(k::BroadcastingFieldOpMap,g::AbstractArray,f::GBlock)
  evaluate(k,g,f)
end

function Fields.return_cache(
  k::BroadcastingFieldOpMap,g::AbstractArray,f::GBlock{A,N}) where {A,N}
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
  GBlock(h,f.touched),l
end

function Fields.evaluate!(cache,k::BroadcastingFieldOpMap,g::AbstractArray,f::GBlock)
  h,l = cache
  @check h.touched == f.touched
  for i in eachindex(f.array)
    if f.touched[i]
      h.array[i] = evaluate!(l[i],k,g,f.array[i])
    end
  end
  h
end

function Fields.return_value(k::BroadcastingFieldOpMap,f::GBlock,g::GBlock)
  evaluate(k,f,g)
end

function Fields.return_cache(k::BroadcastingFieldOpMap,f::GBlock,g::GBlock)
  @notimplemented
end

function Fields.evaluate!(cache,k::BroadcastingFieldOpMap,f::GBlock,g::GBlock)
  @notimplemented
end

function Fields.return_cache(
  k::BroadcastingFieldOpMap,f::GBlock{A,N},g::GBlock{B,N}) where {A,B,N}
  @notimplementedif (k.op != +) && (k.op != -)
  @check size(f) == size(g)
  fi = testvalue(A)
  gi = testvalue(B)
  ci = return_cache(k,fi,gi)
  hi = evaluate!(ci,k,fi,gi)
  a = Array{typeof(hi),N}(undef,size(f.array))
  b = Array{typeof(ci),N}(undef,size(f.array))
  t = map(|,f.touched,g.touched)
  for i in eachindex(f.array)
    if f.touched[i] && g.touched[i]
      b[i] = return_cache(k,f.array[i],g.array[i])
    elseif f.touched[i]
      b[i] = return_cache(BroadcastingFieldOpMap(+),f.array[i])
    elseif g.touched[i]
      b[i] = return_cache(k,g.array[i])
    end
  end
  GBlock(a,t), b
end

function Fields.evaluate!(
  cache,k::BroadcastingFieldOpMap,f::GBlock{A,N},g::GBlock{B,N}) where {A,B,N}
  a,b = cache
  @notimplementedif (k.op != +) && (k.op != -)
  @check size(f) == size(g)
  @check size(a) == size(g)
  for i in eachindex(f.array)
    if f.touched[i] && g.touched[i]
      a.array[i] = evaluate!(b[i],k,f.array[i],g.array[i])
    elseif f.touched[i]
      a.array[i] = evaluate!(b[i],BroadcastingFieldOpMap(+),f.array[i])
    elseif g.touched[i]
      a.array[i] = evaluate!(b[i],k,g.array[i])
    end
  end
  a
end

