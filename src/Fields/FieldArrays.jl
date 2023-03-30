
# Make arrays of field behave like Maps

function return_cache(f::AbstractArray{T},x::Point) where T<:Field
  S = return_type(testitem(f),x)
  cr = CachedArray(zeros(S,size(f)))
  if isconcretetype(T)
    cf = return_cache(testitem(f),x)
  else
    cf = nothing
  end
  cr, cf
end

"""
Implementation of `return_cache` for a array of `Field`.

If the field vector has length `nf` and it is evaluated in one point, it
returns an `nf` vector with the result. If the same array is applied to a
vector of `np` points, it returns a matrix `np` x `nf`.
"""
function evaluate!(c,f::AbstractArray{T},x::Point) where T<:Field
  cr, cf = c
  setsize!(cr,size(f))
  r = cr.array
  if isconcretetype(T)
    for j in eachindex(f)
      @inbounds r[j] = evaluate!(cf,f[j],x)
    end
  else
    for j in eachindex(f)
      @inbounds r[j] = evaluate(f[j],x)
    end
  end
  r
end

function return_cache(f::AbstractArray{T},x::AbstractArray{<:Point}) where T<:Field
  S = return_type(testitem(f),testitem(x))
  cr = CachedArray(zeros(S,(size(x)...,size(f)...)))
  if isconcretetype(T)
    cf = return_cache(f,testitem(x))
  else
    cf = nothing
  end
  cr, cf
end

function evaluate!(c,f::AbstractArray{T},x::AbstractArray{<:Point}) where T<:Field
  cr, cf = c
  setsize!(cr,(size(x)...,size(f)...))
  r = cr.array
  if isconcretetype(T)
    for i in eachindex(x)
      fxi = evaluate!(cf,f,x[i])
      for j in CartesianIndices(f)
        r[i,j] = fxi[j]
      end
    end
  else
    for j in eachindex(f)
      for i in eachindex(x)
        r[i,j] = evaluate(f[j],x[i])
      end
    end
  end
  r
end

function testargs(f::AbstractArray{T},x::Point) where T<:Field
  testargs(testitem(f),x)
end

function testargs(f::AbstractArray{T},x::AbstractArray{<:Point}) where T<:Field
  testargs(testitem(f),x)
end

function test_field_array(f::AbstractArray{<:Field}, x, v, cmp=(==); grad=nothing, gradgrad=nothing)
  test_map(v,f,x;cmp=cmp)
  if grad != nothing
    test_map(grad,Broadcasting(∇)(f),x;cmp=cmp)
  end
  if gradgrad != nothing
    test_map(gradgrad,Broadcasting(∇∇)(f),x;cmp=cmp)
  end
  true
end

# Opening the door to optimize arrays of field gradients

"""
A wrapper that represents the broadcast of `gradient` over an array of fields.
Ng is the number of times the gradient is applied
"""
struct FieldGradientArray{Ng,A,T,N} <: AbstractArray{T,N}
  fa::A
  function FieldGradientArray{Ng}(f::AbstractArray{<:Field}) where Ng
    T = typeof(gradient(testitem(f),Val(Ng)))
    N = ndims(f)
    A = typeof(f)
    new{Ng,A,T,N}(f)
  end
end

function return_value(k::Broadcasting{typeof(∇)},a::AbstractArray{<:Field})
  evaluate(k,a)
end

function return_value(k::Broadcasting{typeof(∇∇)},a::AbstractArray{<:Field})
  evaluate(k,a)
end

function evaluate!(cache,k::Broadcasting{typeof(∇)},a::AbstractArray{<:Field})
  FieldGradientArray{1}(a)
end

function evaluate!(cache,k::Broadcasting{typeof(∇)},a::FieldGradientArray{N}) where N
  FieldGradientArray{N+1}(a.fa)
end

function evaluate!(cache,k::Broadcasting{typeof(∇∇)},a::AbstractArray{<:Field})
  FieldGradientArray{2}(a)
end

function gradient(a::AbstractArray{<:Field})
  msg =
  """\n
  Function gradient (aka ∇) is not defined for arrays of Field objects.
  Use Broadcasting(∇) instead.
  """
  @unreachable msg
end

function ∇∇(a::AbstractArray{<:Field})
  msg =
  """\n
  Double gradient application (aka ∇∇) is not defined for arrays of Field objects.
  Use Broadcasting(∇∇) instead.
  """
  @unreachable msg
end

Base.size(a::FieldGradientArray) = size(a.fa)
Base.axes(a::FieldGradientArray) = axes(a.fa)
Base.getindex(a::FieldGradientArray{Ng},i::Integer) where Ng = gradient(a.fa[i],Val(Ng))
Base.getindex(
   a::FieldGradientArray{Ng,A,T,N},i::Vararg{Integer,N}) where {Ng,A,T,N} = gradient(a.fa[i...],Val(Ng))
Base.IndexStyle(::Type{<:FieldGradientArray{Ng,A}}) where {Ng,A} = IndexStyle(A)

# Optimizing linear_combination.

function linear_combination(a::AbstractVector{<:Number},b::AbstractVector{<:Field})
  column = 1
  LinearCombinationField(a,b,column)
end

struct LinearCombinationField{V,F} <: Field
  values::V
  fields::F
  column::Int
end

for T in (:(Point),:(AbstractVector{<:Point}))
  @eval begin

  function return_cache(a::LinearCombinationField,x::$T)
    cf = return_cache(a.fields,x)
    fx = return_value(a.fields,x)
    v = a.values
    k = LinearCombinationMap(a.column)
    ck = return_cache(k,v,fx)
    cf, ck
  end

  function evaluate!(cache,a::LinearCombinationField,x::$T)
    cf, ck = cache
    fx = evaluate!(cf,a.fields,x)
    v = a.values
    k = LinearCombinationMap(a.column)
    evaluate!(ck,k,v,fx)
  end

  end
end

for op in (:∇,:∇∇)
  @eval begin
    function $op(a::LinearCombinationField)
      fields = Broadcasting($op)(a.fields)
      LinearCombinationField(a.values,fields,a.column)
    end
  end
end

function linear_combination(a::AbstractMatrix{<:Number},b::AbstractVector{<:Field})
  #[ LinearCombinationField(a,b,i) for i in 1:size(a,2) ]
  LinearCombinationFieldVector(a,b)
end

struct LinearCombinationFieldVector{V,F} <: AbstractVector{LinearCombinationField{V,F}}
  values::V
  fields::F
  function LinearCombinationFieldVector(values::AbstractMatrix{<:Number},fields::AbstractVector{<:Field})
    @check size(values,1) == length(fields) """\n
    Incompatible sizes for performing the linear combination

        linear_combination(values,fields) = transpose(values)*fields

    size(values,1) != length(fields)
    """
    V = typeof(values)
    F = typeof(fields)
    new{V,F}(values,fields)
  end
end

Base.size(a::LinearCombinationFieldVector) = (size(a.values,2),)
Base.getindex(a::LinearCombinationFieldVector,i::Integer) = LinearCombinationField(a.values,a.fields,i)
Base.IndexStyle(::Type{<:LinearCombinationField}) = IndexLinear()

for T in (:(Point),:(AbstractVector{<:Point}))
  @eval begin

    function return_cache(a::LinearCombinationFieldVector,x::$T)
      cf = return_cache(a.fields,x)
      fx = return_value(a.fields,x)
      v = a.values
      k = LinearCombinationMap(:)
      ck = return_cache(k,v,fx)
      cf, ck
    end

    function evaluate!(cache,a::LinearCombinationFieldVector,x::$T)
      cf, ck = cache
      fx = evaluate!(cf,a.fields,x)
      v = a.values
      k = LinearCombinationMap(:)
      evaluate!(ck,k,v,fx)
    end

  end
end

for op in (:∇,:∇∇)
  @eval begin
    function evaluate!(cache,k::Broadcasting{typeof($op)},a::LinearCombinationFieldVector)
      fields = Broadcasting($op)(a.fields)
      LinearCombinationFieldVector(a.values,fields)
    end
  end
end

function get_children(n::TreeNode, a::LinearCombinationFieldVector)
  (similar_tree_node(n,a.values),similar_tree_node(n,a.fields))
end

# This is the map that acts on values
struct LinearCombinationMap{T} <: Map
  column::T
  LinearCombinationMap(column::Integer)  = new{typeof(column)}(column)
  LinearCombinationMap(column::Colon)  = new{typeof(column)}(column)
end

function evaluate!(cache,k::LinearCombinationMap{<:Integer},v::AbstractArray,fx::AbstractVector)
  z = zero(return_type(outer,testitem(fx),testitem(v)))
  @check length(fx) == size(v,1)
  @inbounds for i in eachindex(fx)
    # We need to do the product in this way
    # so that the gradient also works
    z += outer(fx[i],v[i,k.column])
  end
  z
end

function return_value(k::LinearCombinationMap{<:Integer},v::AbstractArray,fx::AbstractMatrix)
  if size(fx,2) == size(v,1)
    evaluate(k,v,fx)
  else
    c = return_cache(k,v,fx)
    c.array
  end
end

function return_value(k::LinearCombinationMap{<:Integer},v::AbstractVector,fx::AbstractVector)
  Ta = eltype(v)
  Tb = eltype(fx)
  za = zero(Ta)
  zb = zero(Tb)
  zero( zb⊗za + zb⊗za )
end

function return_cache(k::LinearCombinationMap{<:Integer},v::AbstractArray,fx::AbstractMatrix)
  vf = testitem(fx)
  vv = testitem(v)
  T = typeof( vf⊗vv + vf⊗vv )
  r = zeros(T,size(fx,1))
  CachedArray(r)
end

function evaluate!(cache,k::LinearCombinationMap{<:Integer},v::AbstractArray,fx::AbstractMatrix)
  @check size(fx,2) == size(v,1)
  setsize!(cache,(size(fx,1),))
  r = cache.array
  z = zero(eltype(r))
  @inbounds for p in 1:size(fx,1)
    rp = z
    for i in 1:size(fx,2)
      rp += outer(fx[p,i],v[i,k.column])
    end
    r[p] = rp
  end
  r
end

function evaluate!(cache,k::LinearCombinationMap{Colon},v::AbstractVector,fx::AbstractVector)
  evaluate!(cache,LinearCombinationMap(1),v,fx)
end

function return_value(k::LinearCombinationMap{Colon},v::AbstractVector,fx::AbstractMatrix)
  return_value(LinearCombinationMap(1),v,fx)
end

function return_value(k::LinearCombinationMap{Colon},v::AbstractVector,fx::AbstractVector)
  return_value(LinearCombinationMap(1),v,fx)
end

function return_cache(k::LinearCombinationMap{Colon},v::AbstractVector,fx::AbstractMatrix)
  return_cache(LinearCombinationMap(1),v,fx)
end

function evaluate!(cache,k::LinearCombinationMap{Colon},v::AbstractVector,fx::AbstractMatrix)
  evaluate!(cache,LinearCombinationMap(1),v,fx)
end

function return_cache(k::LinearCombinationMap{Colon},v::AbstractMatrix,fx::AbstractVector)
  vf = testitem(fx)
  vv = testitem(v)
  T = typeof( vf⊗vv + vf⊗vv )
  r = zeros(T,size(v,2))
  CachedArray(r)
end

function evaluate!(cache,k::LinearCombinationMap{Colon},v::AbstractMatrix,fx::AbstractVector)
  @check length(fx) == size(v,1)
  setsize!(cache,(size(v,2),))
  r = cache.array
  @inbounds for j in eachindex(r)
    rj = zero(eltype(r))
    for i in eachindex(fx)
      rj += outer(fx[i],v[i,j])
    end
    r[j] = rj
  end
  r
end

function return_cache(k::LinearCombinationMap{Colon},v::AbstractMatrix,fx::AbstractMatrix)
  vf = testitem(fx)
  vv = testitem(v)
  T = typeof( vf⊗vv + vf⊗vv )
  r = zeros(T,size(fx,1),size(v,2))
  CachedArray(r)
end

function evaluate!(cache,k::LinearCombinationMap{Colon},v::AbstractMatrix,fx::AbstractMatrix)
  @check size(fx,2) == size(v,1)
  setsize!(cache,(size(fx,1),size(v,2)))
  r = cache.array
  @inbounds for p in 1:size(fx,1)
    for j in 1:size(r,2)
      rj = zero(eltype(r))
      for i in 1:size(fx,2)
        rj += outer(fx[p,i],v[i,j])
      end
      r[p,j] = rj
    end
  end
  r
end

# Optimizing transpose
testitem(a::Transpose{<:Field}) = testitem(a.parent)
evaluate!(cache,k::Broadcasting{typeof(∇)},a::Transpose{<:Field}) = transpose(k(a.parent))
evaluate!(cache,k::Broadcasting{typeof(∇∇)},a::Transpose{<:Field}) = transpose(k(a.parent))

return_cache(k::Transpose{<:Field},x::Point) = return_cache(k.parent,x)
evaluate!(cache,k::Transpose{<:Field},x::Point) = transpose(evaluate!(cache,k.parent,x))

return_cache(k::Transpose{<:Field},x::AbstractVector{<:Point}) = return_cache(k.parent,x)
function evaluate!(cache,k::Transpose{<:Field},x::AbstractVector{<:Point})
  TransposeFieldIndices(evaluate!(cache,k.parent,x))
end

struct TransposeMap <: Map end
evaluate!(cache,k::TransposeMap,a::AbstractVector) = transpose(a)
evaluate!(cache,k::TransposeMap,a::AbstractMatrix) = TransposeFieldIndices(a)

"""
Given a matrix `np` x `nf1` x `nf2` result of the evaluation of a field vector
on a vector of points, it returns an array in which the field axes (second and
third axes) are permuted. It is equivalent as `Base.permutedims(A,(1,3,2))`
but more performant, since it does not involve allocations.
"""
struct TransposeFieldIndices{A,T} <: AbstractArray{T,3}
  matrix::A
  function TransposeFieldIndices(matrix::AbstractMatrix{T}) where T
    A = typeof(matrix)
    new{A,T}(matrix)
  end
end

function TransposeFieldIndices{A,T}(::UndefInitializer,shape::NTuple{3,Integer}) where {A,T}
  TransposeFieldIndices(similar(A,(shape[1],shape[3])))
end
Base.size(a::TransposeFieldIndices) = (size(a.matrix,1),1,size(a.matrix,2))
Base.axes(a::TransposeFieldIndices) = (axes(a.matrix,1),Base.OneTo(1),axes(a.matrix,2))
Base.IndexStyle(::Type{<:TransposeFieldIndices{A}}) where A = IndexStyle(A)
Base.getindex(a::TransposeFieldIndices,i::Integer,j::Integer,k::Integer) = a.matrix[i,k]
Base.getindex(a::TransposeFieldIndices,i::Integer) = a.matrix[i]
Base.setindex!(a::TransposeFieldIndices,v,i::Integer,j::Integer,k::Integer) = (a.matrix[i,k] = v)
Base.setindex!(a::TransposeFieldIndices,v,i::Integer) = (a.matrix[i] = v)


# Integration

"""
Integration of a given array of fields in the "physical" space
"""
function integrate(a::AbstractArray{<:Field},x::AbstractVector{<:Point},w::AbstractVector{<:Real})
  cache = return_cache(integrate,a,x,w)
  evaluate!(cache,integrate,a,x,w)
end

"""
Integration of a given array of fields in the "reference" space
"""
function integrate(a::AbstractArray{<:Field},q::AbstractVector{<:Point},w::AbstractVector{<:Real},j::Field)
  cache = return_cache(integrate,a,q,w,j)
  evaluate!(cache,integrate,a,q,w,j)
end

# Broadcast operations

function return_value(k::Broadcasting{<:Operation},args::Union{Field,AbstractArray{<:Field}}...)
  BroadcastOpFieldArray(k.f.op,args...)
end

function evaluate!(cache,k::Broadcasting{<:Operation},args::Union{Field,AbstractArray{<:Field}}...)
  BroadcastOpFieldArray(k.f.op,args...)
end

"""
Type that represents a broadcast operation over a set of `AbstractArray{<:Field}`.
The result is a sub-type of `AbstractArray{<:Field}`
"""
struct BroadcastOpFieldArray{O,T,N,A} <: AbstractArray{T,N}
  op::O
  args::A
  function BroadcastOpFieldArray(op,args::Union{Field,AbstractArray{<:Field}}...)
    fs = map(testitem,args)
    T = return_type(Operation(op),fs...)
    s = map(size,args)
    bs = Base.Broadcast.broadcast_shape(s...)
    N = length(bs)
    A = typeof(args)
    O = typeof(op)
    new{O,T,N,A}(op,args)
  end
end

Base.size(a::BroadcastOpFieldArray) = Base.Broadcast.broadcast_shape(map(size,a.args)...)
Base.axes(a::BroadcastOpFieldArray) = Base.Broadcast.broadcast_shape(map(axes,a.args)...)
Base.IndexStyle(::Type{<:BroadcastOpFieldArray}) = IndexLinear()
Base.getindex(a::BroadcastOpFieldArray,i::Integer) = broadcast(Operation(a.op),a.args...)[i]
function testitem(a::BroadcastOpFieldArray)
  fs = map(testitem,a.args)
  return_value(Operation(a.op),fs...)
end

for T in (:(Point),:(AbstractArray{<:Point}))
  @eval begin

    function return_cache(f::BroadcastOpFieldArray,x::$T)
      cfs = map(fi -> return_cache(fi,x),f.args)
      rs = map(fi -> return_value(fi,x),f.args)
      bm = BroadcastingFieldOpMap(f.op)
      r = return_cache(bm,rs...)
      r, cfs
    end

    function evaluate!(c,f::BroadcastOpFieldArray,x::$T)
      r, cfs = c
      rs = map((ci,fi) -> evaluate!(ci,fi,x),cfs,f.args)
      bm = BroadcastingFieldOpMap(f.op)
      evaluate!(r,bm,rs...)
    end

  end
end

# With this type we mark that we are doing Broadcasting(op) on the result of evaluating Fields/FieldArrays
# This allow us to do some optimizations for block arrays that are only true in this context, not in a
# general Broadcasting operation.
struct BroadcastingFieldOpMap{F} <: Map
  op::F
end

return_value(a::BroadcastingFieldOpMap,args...) = return_value(Broadcasting(a.op),args...)
return_cache(a::BroadcastingFieldOpMap,args...) = return_cache(Broadcasting(a.op),args...)
evaluate!(cache,a::BroadcastingFieldOpMap,args...) = evaluate!(cache,Broadcasting(a.op),args...)

return_value(a::BroadcastingFieldOpMap,args::AbstractArray...) = return_value(Broadcasting(a.op),args...)
return_cache(a::BroadcastingFieldOpMap,args::AbstractArray...) = return_cache(Broadcasting(a.op),args...)
evaluate!(cache,a::BroadcastingFieldOpMap,args::AbstractArray...) = evaluate!(cache,Broadcasting(a.op),args...)

# Follow optimizations are very important to achieve performance

function evaluate!(
  cache,
  f::BroadcastingFieldOpMap,
  a::AbstractArray{T,N},
  b::AbstractArray{S,N}) where {T,S,N}

  @check size(a) == size(b) || (length(a)==0 && length(b)==0)
  setsize!(cache,size(a))
  r = cache.array
  for i in eachindex(a)
    r[i] = f.op(a[i],b[i])
  end
  r
end

function evaluate!(
  cache,
  f::BroadcastingFieldOpMap,
  a::AbstractMatrix,
  b::AbstractArray{S,3} where S)

  @check size(a,1) == size(b,1)
  @check size(b,2) == 1 || size(b,1) == 0
  np, ni = size(a)
  nj = size(b,3)
  setsize!(cache,(np,ni,nj))
  r = cache.array
  for j in 1:nj
    for p in 1:np
      bpj = b[p,1,j]
      for i in 1:ni
        r[p,i,j] = f.op(a[p,i],bpj)
      end
    end
  end
  r
end

function evaluate!(
  cache,
  f::BroadcastingFieldOpMap,
  b::AbstractArray{S,3} where S,
  a::AbstractMatrix)

  @check size(a,1) == size(b,1)
  @check size(b,2) == 1 || size(b,1) == 0
  np, ni = size(a)
  nj = size(b,3)
  setsize!(cache,(np,ni,nj))
  r = cache.array
  for p in 1:np
    for j in 1:nj
      bpj = b[p,1,j]
      for i in 1:ni
        r[p,i,j] = f.op(bpj,a[p,i])
      end
    end
  end
  r
end

function evaluate!(
  cache,
  f::BroadcastingFieldOpMap,
  a::AbstractVector,
  b::AbstractMatrix)

  @check size(a,1) == size(b,1)
  np, ni = size(b)
  setsize!(cache,(np,ni))
  r = cache.array
  for p in 1:np
    ap = a[p]
    for i in 1:ni
      r[p,i] = f.op(ap,b[p,i])
    end
  end
  r
end

function evaluate!(
  cache,
  f::BroadcastingFieldOpMap,
  b::AbstractMatrix,
  a::AbstractVector)

  @check size(a,1) == size(b,1)
  np, ni = size(b)
  setsize!(cache,(np,ni))
  r = cache.array
  for p in 1:np
    ap = a[p]
    for i in 1:ni
      r[p,i] = f.op(b[p,i],ap)
    end
  end
  r
end

function evaluate!(
  cache,
  f::BroadcastingFieldOpMap,
  a::AbstractVector,
  b::AbstractArray{S,3} where S)

  @check size(a,1) == size(b,1)
  np, ni, nj = size(b)
  setsize!(cache,(np,ni,nj))
  r = cache.array
  for p in 1:np
    ap = a[p]
    for j in 1:nj
      for i in 1:ni
        r[p,i,j] = f.op(ap,b[p,i,j])
      end
    end
  end
  r
end

function evaluate!(
  cache,
  f::BroadcastingFieldOpMap,
  b::AbstractArray{S,3} where S,
  a::AbstractVector)

  @check size(a,1) == size(b,1)
  np, ni, nj = size(b)
  setsize!(cache,(np,ni,nj))
  r = cache.array
  for p in 1:np
    ap = a[p]
    for j in 1:nj
      for i in 1:ni
        r[p,i,j] = f.op(b[p,i,j],ap)
      end
    end
  end
  r
end

# Gradient of the sum
for op in (:+,:-)
  @eval begin
    function evaluate!(cache,::Broadcasting{typeof(∇)},a::BroadcastOpFieldArray{typeof($op)})
      f = a.args
      g = map( Broadcasting(∇), f)
      Broadcasting(Operation($op))(g...)
    end
  end
end

# Gradient of the product
for op in (:*,:⋅,:⊙,:⊗)
  @eval begin
    function evaluate!(cache,::Broadcasting{typeof(∇)},a::BroadcastOpFieldArray{typeof($op)})
      f = a.args
      @notimplementedif length(f) != 2
      f1, f2 = f
      g1, g2 = map(Broadcasting(∇), f)
      k(F1,F2,G1,G2) = product_rule($op,F1,F2,G1,G2)
      Broadcasting(Operation(k))(f1,f2,g1,g2)
    end
  end
end
