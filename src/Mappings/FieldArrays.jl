# Arrays of Fields

# Default implementation of evaluate! for arrays of fields

function return_cache(f::AbstractArray{<:Field},x::Point)
  cfi = return_cache(first(f),x)
  T = return_type(first(f),x)
  sf = size(f)
  _cf = zeros(T,sf)
  cf = CachedArray(_cf)
  cf, cfi
end

function evaluate!(c,f::AbstractArray{<:Field},x::Point)
  cf, cfi = c
  sf = size(f)
  setsize!(cf,sf)
  for i in eachindex(f)
    @inbounds cf.array[i] = evaluate!(cfi,f[i],x)
  end
  cf.array
end

function return_cache(f::AbstractArray{<:BroadcastField},x::AbstractArray{<:Point})
  fi = first(f).field
  xi = first(x)
  cfi = return_cache(fi,xi)
  T = return_type(fi,xi)
  sf = size(f)
  sx = size(x)
  _cf = zeros(T,sx...,sf...)
  cf = CachedArray(_cf)
  cf, cfi
end

function evaluate!(c,f::AbstractArray{<:BroadcastField},x::AbstractArray{<:Point})
  cf, cfi = c
  sf = size(f)
  sx = size(x)
  s = (sx...,sf...)
  setsize!(cf,s)
  for j in CartesianIndices(f)
    for i in CartesianIndices(x)
      @inbounds cf.array[i,j] = evaluate!(cfi,f[j].field,x[i])
    end
  end
  cf.array
end

# Gradient

# We can create optimised versions for concrete types
gradient(fa::AbstractArray{<:Field}) = gradient.(fa)

# Operations

# @santiagobadia : Just to check, think further
function *(A::AbstractMatrix{<:Field}, B::AbstractMatrix{<:Field})
  # TS = promote_op(matprod, eltype(A), eltype(B))
  TS = Field
  mul!(similar(B, TS, (size(A,1), size(B,2))), A, B)
end

#Compute operations with arrays efficiently

evaluate!(cache,op::Operation,x::AbstractArray{<:Field}...) = GenericField.(OperationMapping.(op.op,x))

struct OperationArray{T,N,S<:Field} <: AbstractArray{T,N}
  op::T
  res::AbstractArray{S,N}
  args
  function OperationArray(op::Operation,args::Union{Field,AbstractArray{<:Field}}...)
    res = op.op(args...)
    S = eltype(res)
    new{typeof(op),ndims(res),S}(op,res,args)
  end
end

function return_cache(oa::OperationArray,x)
  ci = return_caches(oa.args,x)
  sr = size(oa.res)
  # T = return_type(oa.op.op,return_types(oa.args,x)...)
  # cr = zero(T(undef,sr))
  ri = map(a->testitem(a,x),oa.args)
  r = testitem(oa.op.op,ri...)
  r, ci
  # Missing Cached arrays r cr
end

function evaluate!(c,f::OperationArray,x)
  r, ci = c
  # @santiagobadia : Not sure where the allocation is ???
  # 3 allocations independent of array size, number operands,
  # types of operands, etc
  _inplace!(f.op.op,r,_evaluate!(ci,x,f.args)...)
  r
end

function _evaluate!(c,x,fs)
  @inbounds for i in 1:length(fs)
    evaluate!(c[i],fs[i],x)
  end
  map(ci->ci[1].array,c)
end

Base.size(oa::OperationArray) = Base.size(oa.res)
Base.axes(oa::OperationArray) = Base.axes(oa.res)
Base.IndexStyle(oa::OperationArray) where A = Base.IndexStyle(oa.res)
Base.eltype(oa::OperationArray) = Base.eltype(oa.res)
Base.iterate(oa::OperationArray) = Base.iterate(oa.res)
Base.getindex(oa::OperationArray,i) = Base.getindex(oa.res,i)

# @santiagobadia : I have created these inplace array operations

@inline function _inplace!(op::Union{typeof(+),typeof(-)},r::AbstractArray,a::AbstractArray)
  for i in eachindex(r)
    r[i] = op(a[i])
  end
end

@inline function _inplace!(op::Union{typeof(+),typeof(-)},r::AbstractArray,a::AbstractArray,b::AbstractArray)
  for i in eachindex(r)
    r[i] = op(a[i],b[i])
  end
end

@inline function _inplace!(op::Union{typeof(+),typeof(-)},r::AbstractArray,a::AbstractArray,b::AbstractArray,c::AbstractArray...)
  _inplace!(op,r,a,b)
  _inplace!(op,r,r,c...)
end

@inline function _inplace!(::typeof(*),r::AbstractArray,a::AbstractArray) end

@inline function _inplace!(::typeof(*),r::AbstractArray,a::AbstractArray,b::AbstractArray)
  mul!(r,a,b)
end

@inline function _inplace!(::typeof(*),r::AbstractArray,a::AbstractArray,b::AbstractArray,c::AbstractArray...)
  @notimplemented
end

# Testers

function test_field_array(f,p,nf;grad=false,hessian=false)
  fa = fill(f,nf...)
  fp = evaluate(f,p)
  fap = fill(fp,nf...)
  test_field(fa,p,fap)
  if hessian
    test_field_array(∇(f),p,nf;grad=true)
  elseif grad
    test_field_array(∇(f),p,nf)
  end
  fa, p
end

function test_broadcast_field_array(f,p,nf,np;grad=false,hessian=false)
  bf = BroadcastField(f)
  bfa = fill(bf,nf...)
  fp = evaluate(f,p)
  x = fill(p,np...)
  fax = fill(fp,np...,nf...)
  test_field(bfa,x,fax)
  if hessian
    test_broadcast_field_array(∇(f),p,nf,np;grad=true)
  elseif grad
    test_broadcast_field_array(∇(f),p,nf,np)
  end
  bfa, x
end
