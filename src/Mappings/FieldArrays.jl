# Arrays of Fields

# Default implementation of evaluate! for arrays of fields

function return_cache(f::AbstractArray{<:Field},x::Point)
  # @santiagobadia : Be careful, we could have different fields
  _cfs = map(fi -> return_cache(fi,x),f)
  cfs = _cfs # Not sure about CachedArray here, it does not have much sense
  # @santiagobadia : We should check all return the same
  # Anyway, this can be overwritten by concrete AbstractArray{<:Field}
  T = return_type(first(f),x)
  sf = size(f)
  _cf = zeros(T,sf)
  cf = CachedArray(_cf)
  cf, cfs
end

function evaluate!(c,f::AbstractArray{<:Field},x::Point)
  cf, cfs = c
  sf = size(f)
  # Expect error here if f changes size and not type-stable
  setsize!(cf,sf)
  for i in eachindex(f)
    # @ santiagobadia : Not allocation free
    @inbounds cf.array[i] = evaluate!(cfs[i],f[i],x)
  end
  cf.array
end

function return_cache(f::AbstractArray{<:Field},x::AbstractArray{<:Point})
  fi = first(f)
  xi = first(x)
  # @santiagobadia : Same comment as above
  cfi = return_cache(fi,xi)
  T = return_type(fi,xi)
  sf = size(f)
  sx = size(x)
  _cf = zeros(T,sx...,sf...)
  cf = CachedArray(_cf)
  cf, cfi
end

function evaluate!(c,f::AbstractArray{<:Field},x::AbstractArray{<:Point})
  cf, cfi = c
  sf = size(f)
  sx = size(x)
  s = (sx...,sf...)
  setsize!(cf,s)
  for j in CartesianIndices(f)
    for i in CartesianIndices(x)
      @inbounds cf.array[i,j] = evaluate!(cfi,f[j],x[i])
    end
  end
  cf.array
end

# Gradient

# We can create optimised versions for concrete types (but this is very fast)

gradient(fa::AbstractArray{<:Field}) = gradient.(fa)

# Operations

# @santiagobadia : Probably improve it in the future, with some promote method
# for Fields, not really to be used in practise

const FieldMatrix{T,S} = Union{Transpose{T,S},S} where {T<:Field,S<:AbstractMatrix{<:Field}}

function (*)(A::FieldMatrix, B::FieldMatrix)
  TS = Field
  mul!(similar(B, TS, (size(A,1), size(B,2))), A, B)
end

# Transpose Field Array

return_cache(f::Transpose{<:Field},x::Point) = return_cache(transpose(f),x)
return_cache(f::Transpose{<:Field},x::AbstractArray{<:Point}) = return_cache(transpose(f),x)
evaluate!(c,f::Transpose{<:Field},x::Point) = transpose(evaluate!(c,transpose(f),x))

function evaluate!(c,f::Transpose{<:Field},x::AbstractArray{<:Point})
  _r = evaluate!(c,transpose(f),x)
  # @santiagobadia : Other ideas?
  r = permutedims(reshape(_r,size(_r)...,1),[1,3,2])
end

# Compute operations with arrays efficiently

struct OperationArray{T,N,S<:Field} <: AbstractArray{S,N}
  op::T
  res::AbstractArray{S,N}
  args
  function OperationArray(op::Operation,args::Union{Field,AbstractArray{<:Field}}...)
    res = op.op(args...)
    S = eltype(res)
    new{typeof(op.op),ndims(res),S}(op.op,res,args)
  end
end

# @santiagobadia : We can change this, but not much important
# Probably reshape it eliminating the first dimension
function return_cache(f::OperationArray,x::Point)
  return_cache(f,[x])
end

function return_cache(f::OperationArray{typeof(transpose)},x::AbstractArray{<:Point})
  @assert length(f.args) == 1
  return_cache(f.args[1],x)
end

# @santiagobadia : Is this what we want here?
function evaluate!(c,f::OperationArray{typeof(transpose)},x::AbstractArray{<:Point})
  transpose(evaluate!(c,f.args...,x))
end

# @santiagobadia : Do it with eval for + -
# for op in (:+,:-,:*,:⋅,:inv,:det,:transpose)
# eval ...
function return_cache(f::OperationArray{typeof(+)},x::AbstractArray{<:Point})
  cfs = return_caches(f.args,x)
  rs = map(fi -> evaluate(fi,first(x)),f.args)
  cr = CachedArray(f.op(rs...))
  cr, cfs
end

# @santiagobadia :  How can we implement it allocation free?
# How can we compute new dim if nfields is changing
function evaluate!(c,f::OperationArray{typeof(+)},x::AbstractArray{<:Point})
  cr, cfs = c
  ri = map((ci,fi)->evaluate!(ci,fi,x),cfs,f.args)
  f.op(ri...)
end

function return_cache(f::OperationArray{typeof(-)},x::AbstractArray{<:Point})
  rs = return_caches(f.args,x)
end

function evaluate!(c,f::OperationArray{typeof(-)},x::AbstractArray{<:Point})
  ri = map((ci,fi)->evaluate!(ci,fi,x),c,f.args)
  f.op(ri...)
end

Base.ndims(::Field) = 0

function (op::Operation{typeof(*)})(A::AbstractArray{<:Field},B::AbstractArray{<:Field})
  r = A*B
  if ndims(r) == 0
    return DotArray(A,B)
  else
    return OperationArray(op,A,B)
  end
end

(op::Operation{typeof(⋅)})(A::AbstractArray{<:Field},B::AbstractArray{<:Field}) = DotArray(A,B)

# Result of array operation is a scalar

struct DotArray <: Field
  f
  g
end

# inner_product(f,g) = DotArray(f,g)

function return_cache(f::DotArray,x::AbstractArray{<:Point})
  f1 = f.f
  f2 = f.g
  c1 = return_cache(f1,x)
  c2 = return_cache(f2,x)
  r1 = evaluate(f1,first(x))
  r2 = evaluate(f2,first(x))
  r = CachedArray(zeros(typeof(r1⋅r2),length(x)))
  r, (c1, c2)
end

function evaluate!(c,f::DotArray,x::AbstractArray{<:Point})
  cr, cfs = c
  c1, c2 = cfs
  f1 = f.f
  f2 = f.g
  setsize!(cr,size(x))
  r = cr.array
  r1 = evaluate!(c1,f1,x)
  r2 = evaluate!(c2,f2,x)
  for i in 1:length(x)
    r[i] = selectdim(r1,1,i)⋅selectdim(r2,1,i)
  end
  r
end

# Result of array operation is an array

# matmul(A,B) = OperationArray(Operation(*),A,B)

function return_cache(f::OperationArray{typeof(*)},x::Point)
  @assert length(f.args) == 2
  f1, f2 = f.args
  c1 = return_cache(f1,x)
  c2 = return_cache(f2,x)
  r1 = evaluate!(c1,f1,first(x))
  r2 = evaluate!(c2,f2,first(x))
  r = CachedArray(zeros(eltype(r1*r2),length(f1),length(f1)))
  r, (c1,c2)
end

function return_cache(f::OperationArray{typeof(*)},x::AbstractArray{<:Point})
  @assert length(f.args) == 2
  f1, f2 = f.args
  c1 = return_cache(f1,first(x))
  c2 = return_cache(f2,first(x))
  r1 = evaluate!(c1,f1,first(x))
  r2 = evaluate!(c2,f2,first(x))
  _r = r1*r2
  r = CachedArray(zeros(eltype(_r),length(x),size(_r)...))
  r, (c1,c2)
end

function evaluate!(c,f::OperationArray{typeof(*)},x::Point)
  cr, cfs = c
  r = cr.array
  c1, c2 = cfs
  f1, f2 = f.args
  r1 = evaluate!(c1,f1,x)
  r2 = evaluate!(c2,f2,x)
  r = r1*r2
end
function evaluate!(c,f::OperationArray{typeof(*)},x::AbstractArray{<:Point})
  cr, cfs = c
  r = cr.array
  c1, c2 = cfs
  f1, f2 = f.args
  for i in eachindex(x)
    xi = x[i]
    r1 = evaluate!(c1,f1,xi)
    r2 = evaluate!(c2,f2,xi)
    _r = r1*r2
    r[i,:,:] = _r
  end
  r
end

# Linear combination

# @santiagobadia : Pobably not really needed...
# linear_combination(f,v) = GenericField.(v)⋅f
linear_combination(f,v) = LinearCombination(f,v)

struct LinearCombination <: Field
  base
  vals
end

function return_cache(lc::LinearCombination,x::AbstractArray{<:Point})
  f = lc.base
  v = lc.vals
  cf = return_cache(f,x)
  rf = evaluate(f,first(x))
  cr = CachedArray(zeros(typeof(rf⋅v),length(x)))
  cr, cf
end

function evaluate!(c,lc::LinearCombination,x::AbstractArray{<:Point})
  cr, cf = c
  f = lc.base
  v = lc.vals
  rf = evaluate!(cf,f,x)
  setsize!(cr,(length(x),))
  for i in eachindex(x)
    cr.array[i] = rf[i,:]⋅v
  end
  cr.array
end

# AbstractArray API for OperationArray

Base.size(oa::OperationArray) = Base.size(oa.res)
Base.axes(oa::OperationArray) = Base.axes(oa.res)
Base.IndexStyle(oa::OperationArray) where A = Base.IndexStyle(oa.res)
Base.eltype(oa::OperationArray) = Base.eltype(oa.res)
Base.iterate(oa::OperationArray) = Base.iterate(oa.res)
Base.getindex(oa::OperationArray,i) = Base.getindex(oa.res,i)


evaluate!(cache,op::Operation,x::AbstractArray{<:Field}...) = OperationArray(op,x...)
evaluate!(cache,op::Operation,x::Union{Field,AbstractArray{<:Field}}...) = OperationArray(op,x...)

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
  test_field_array(f,p,nf,grad=grad,hessian=hessian)
  fa = fill(f,nf...)
  fp = evaluate(f,p)
  x = fill(p,np...)
  fax = fill(fp,np...,nf...)
  test_field(fa,x,fax)
  if hessian
    test_broadcast_field_array(∇(f),p,nf,np;grad=true)
  elseif grad
    test_broadcast_field_array(∇(f),p,nf,np)
  end
  fa, x
end

function test_operation_field_array(op,x,fs...)
  fop = Operation(op)
  fb = op(fs...)
  fba = fop(fs...)
  if fb isa AbstractArray
    @test eltype(fb) <: Field
    @test fb isa AbstractArray{<:Field}
    @test fba isa OperationArray
    @test fba.res == fb
  else
    @test fb isa Field
  end
  @test evaluate(fb,x) == evaluate(fba,x)
  fba
end
