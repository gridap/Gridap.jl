# Default implementation of evaluate! for arrays of fields

@inline function return_cache(f::AbstractArray{<:Field},x)
  _default_return_cache(f,x)
end

@inline function evaluate!(c,f::AbstractArray{<:Field},x)
  _default_evaluate!(c,f,x)
end

function _default_return_cache(f::AbstractArray{T},x::Point) where T<:Field
  S = return_type(first(f),first(x))
  cr = CachedArray(zeros(S,size(f)))
  if isconcretetype(T)
    cf = return_cache(first(f),first(x))
  else
    cf = nothing
  end
  cr, cf
end

function _default_evaluate!(c,f::AbstractArray{T},x::Point) where T<:Field
  cr, cf = c
  setsize!(cr,size(f))
  if isconcretetype(T)
    for j in eachindex(f)
      @inbounds cr.array[j] = evaluate!(cf,f[j],x)
    end
  else
    for j in eachindex(f)
      @inbounds cr.array[j] = evaluate(f[j],x)
    end
  end
  cr.array
end

function _default_return_cache(f::AbstractArray{T},x::AbstractArray{<:Point}) where T<:Field
  S = return_type(first(f),first(x))
  cr = CachedArray(zeros(S,(size(x)...,size(f)...)))
  if isconcretetype(T)
    cf = return_cache(first(f),first(x))
  else
    cf = nothing
  end
  cr, cf
end

function _default_evaluate!(c,f::AbstractArray{T},x::AbstractArray{<:Point}) where T<:Field
  cr, cf = c
  setsize!(cr,(size(x)...,size(f)...))
  if isconcretetype(T)
    for j in eachindex(f)
      for i in eachindex(x)
        @inbounds cr.array[i,j] = evaluate!(cf,f[j],x[i])
      end
    end
  else
    for j in eachindex(f)
      for i in eachindex(x)
        @inbounds cr.array[i,j] = evaluate(f[j],x[i])
      end
    end
  end
  cr.array
end

# Broadcasting operations

# Non-performant default implementation

@inline return_cache(b::BroadcastMapping,args::Union{Field,AbstractArray{<:Field}}...) = nothing

@inline function evaluate!(cache,b::BroadcastMapping,args::Union{Field,AbstractArray{<:Field}}...)
  broadcast(b.f,args...)
end

# @santiagobadia : I am not sure we want this syntatic sugar
# @santiagobadia : Here if we put AbstractArray it goes to the + in Julia for Array...
# for op in (:+,:-,:*,:⋅)
# @eval $op(a::Array{<:Field},b::Array{<:Field}) = BroadcastMapping(Operation($op))(a,b)
# @eval $op(a::AbstractArray{<:Field},b::AbstractArray{<:Field}) = BroadcastMapping(Operation($op))(a,b)
# end

# Transpose

@inline function Base.transpose(a::AbstractVector{<:Field})
  TransposeFieldVector(a)
end

struct TransposeFieldVector{B,T} <: AbstractArray{T,2}
  basis::B
  function TransposeFieldVector(basis)
    @assert ndims(basis) == 1
    T = eltype(basis)
    new{typeof(basis),T}(basis)
  end
end

@inline Base.size(a::TransposeFieldVector) = (1,length(a.basis))
@inline Base.axes(a::TransposeFieldVector) = (Base.OneTo(1),axes(a.basis,1))
@inline Base.getindex(a::TransposeFieldVector,i::Integer,j::Integer) = a.basis[j,i]
@inline Base.getindex(a::TransposeFieldVector,i::Integer) = a.basis[i]
@inline Base.IndexStyle(::TransposeFieldVector{A}) where A = IndexStyle(A)

@inline return_cache(a::TransposeFieldVector,x) = return_cache(a.basis,x)

function evaluate!(cache,a::TransposeFieldVector,x)#::AbstractArray{<:Point})
  M = evaluate!(cache,a.basis,x)
  transpose_field_indices(M)
end

@inline transpose_field_indices(M::AbstractVector) = transpose(M)
@inline transpose_field_indices(M::AbstractMatrix) = TransposeFieldIndices(M)

struct TransposeFieldIndices{A,T} <: AbstractArray{T,3}
  matrix::A
  @inline function TransposeFieldIndices(matrix::AbstractMatrix{T}) where T
    A = typeof(matrix)
    new{A,T}(matrix)
  end
end

@inline Base.size(a::TransposeFieldIndices) = (size(a.matrix,1),1,size(a.matrix,2))
@inline Base.axes(a::TransposeFieldIndices) = (axes(a.matrix,1),Base.OneTo(1),axes(a.matrix,2))
@inline Base.IndexStyle(::Type{<:TransposeFieldIndices{A}}) where A = IndexStyle(A)
@inline Base.getindex(a::TransposeFieldIndices,i::Integer,j::Integer,k::Integer) = a.matrix[i,k]
@inline Base.getindex(a::TransposeFieldIndices,i::Integer) = a.matrix[i]
@inline Base.setindex!(a::TransposeFieldIndices,v,i::Integer,j::Integer,k::Integer) = (a.matrix[i,k] = v)
@inline Base.setindex!(a::TransposeFieldIndices,v,i::Integer) = (a.matrix[i] = v)

# Broadcast operations

@inline (b::BroadcastMapping{<:Operation})(args::Union{Field,AbstractArray{<:Field}}...) = BroadcastOpFieldArray(b.f.op,args...)

struct BroadcastOpFieldArray{O,T,S<:Field,N} <: AbstractArray{S,N}
  op::O
  args::T
  function BroadcastOpFieldArray(op,args::Union{Field,AbstractArray{<:Field}}...)
    fs = map(first,args)
    S = typeof(op(fs...))
    s = map(size,args)
    bs = Base.Broadcast.broadcast_shape(s...)
    N = length(bs)
    T = typeof(args)
    O = typeof(op)
    new{O,T,S,N}(op,args)
  end
end

@inline Base.size(a::BroadcastOpFieldArray) = Base.Broadcast.broadcast_shape(map(size,a.args)...)
@inline Base.IndexStyle(::Type{<:BroadcastOpFieldArray}) = IndexLinear
@inline Base.getindex(a::BroadcastOpFieldArray,I...) = broadcast(a.op,a.args...)[I...]

function return_cache(f::BroadcastOpFieldArray,x)
  cfs = map(fi -> return_cache(fi,x),f.args)
  # cfs = return_caches(f.args,x)
  rs = map((ci,fi) -> evaluate!(ci,fi,x),cfs,f.args)
  bm = BroadcastMapping(f.op)
  r = return_cache(bm,rs...)
  r, cfs
end

function evaluate!(c,f::BroadcastOpFieldArray,x)
  r, cfs = c
  rs = map((ci,fi) -> evaluate!(ci,fi,x),cfs,f.args)
  bm = BroadcastMapping(f.op)
  evaluate!(r,bm,rs...)
  r.array
end

# Dot product vectors

struct DotOpFieldVectors{F,G} <: Field
  f::F
  g::G
end

@inline *(f::TransposeFieldVector,g::AbstractVector{<:Field}) = DotOpFieldVectors(f.basis,g)
@inline ⋅(f::TransposeFieldVector,g::AbstractVector{<:Field}) = DotOpFieldVectors(f.basis,g)
@inline ⋅(g::AbstractVector{<:Field},f::AbstractVector{<:Field}) = DotOpFieldVectors(f,g)
@inline ⋅(g::AbstractVector{<:Number},f::AbstractVector{<:Field}) = DotOpFieldVectors(f,GenericField.(g))
@inline ⋅(f::AbstractVector{<:Field},g::AbstractVector{<:Number}) = DotOpFieldVectors(f,GenericField.(g))

function return_cache(f::DotOpFieldVectors,x::Point)
  f1 = f.f
  f2 = f.g
  c1 = return_cache(f1,x)
  c2 = return_cache(f2,x)
  c1, c2
end

function evaluate!(c,f::DotOpFieldVectors,x::Point)
  cfs = c
  c1, c2 = cfs
  f1 = f.f
  f2 = f.g
  r = cr.array
  r1 = evaluate!(c1,f1,x)
  r2 = evaluate!(c2,f2,x)
  r1⋅r2
end

function return_cache(f::DotOpFieldVectors,x::AbstractArray{<:Point})
  f1 = f.f
  f2 = f.g
  c1 = return_cache(f1,x)
  c2 = return_cache(f2,x)
  r1 = evaluate(f1,first(x))
  r2 = evaluate(f2,first(x))
  r = CachedArray(zeros(typeof(r1⋅r2),length(x)))
  r, (c1, c2)
end

function evaluate!(c,f::DotOpFieldVectors,x::AbstractArray{<:Point})
  cr, cfs = c
  c1, c2 = cfs
  f1 = f.f
  f2 = f.g
  setsize!(cr,size(x))
  r = cr.array
  r1 = evaluate!(c1,f1,x)
  r2 = evaluate!(c2,f2,x)
  for i in 1:length(x)
    @inbounds r[i] = view(r1,i,:)⋅view(r2,i,:)
  end
  r
end

# linear combination

# Function to be used for dispatch in MappedArray
@inline linear_combination(a::AbstractVector{<:Field},b::AbstractArray{<:Number}) = transpose(a)*b
@inline linear_combination(a::AbstractArray{<:Number},b::AbstractVector{<:Field}) = transpose(b)*a

# mat product
#
# The only user API needed
# transpose(i_to_f)*i_to_vals
# transpose(i_to_f)*ij_to_vals
#
@inline function Base.:*(a::TransposeFieldVector,b::AbstractVector{<:Number})
  LinearCombinationField(a.basis,b)
end

@inline function Base.:*(a::TransposeFieldVector,b::AbstractMatrix{<:Number})
  LinearCombinationField(a.basis,b)
end

# These ones are not needed
Base.:*(a::AbstractMatrix{<:Field},b::AbstractVector{<:Number}) = @notimplemented
Base.:*(a::AbstractMatrix{<:Field},b::AbstractMatrix{<:Number}) = @notimplemented

struct LinearCombinationField{B,V} <: Field
  basis::B
  values::V
end

@inline return_cache(a::LinearCombinationField,x::Point) = _lincomb_return_cache(a,x)
@inline evaluate!(cache,a::LinearCombinationField,x::Point) = _lincomb_evaluate!(cache,a,x)

@inline return_cache(a::LinearCombinationField,x::AbstractArray{<:Point}) = _lincomb_return_cache(a,x)
@inline evaluate!(cache,a::LinearCombinationField,x::AbstractArray{<:Point}) = _lincomb_evaluate!(cache,a,x)

function _lincomb_return_cache(a::LinearCombinationField,x)#::AbstractArray{<:Point})
  cb = return_cache(a.basis,x)
  bx = evaluate!(cb, a.basis, x)
  v = a.values
  cr = return_cache(LinCombVal(),bx,v)
  cb, cr
end

function _lincomb_evaluate!(cache,a::LinearCombinationField,x)#::AbstractArray{<:Point})
  cb, cr = cache
  bx = evaluate!(cb,a.basis,x)
  evaluate!(cr,LinCombVal(),bx,a.values)
end

# Composition

struct CompositionFieldArrayField{T,N,A,B} <:AbstractArray{T,N}
  f::A
  g::B
  function CompositionFieldArrayField(f::AbstractArray{<:Field},g::Field)
    s = size(f)
    T = typeof(first(f)∘g)
    N = length(s)
    A = typeof(f)
    B = typeof(g)
    new{T,N,A,B}(f,g)
  end
end

@inline (b::BroadcastMapping{typeof(∘)})(f::AbstractArray{<:Field},g::Field) = CompositionFieldArrayField(f,g)

function return_cache(fa::CompositionFieldArrayField,x)#::AbstractArray{<:Point})
  cg = return_cache(fa.g,x)
  rg = evaluate!(cg,fa.g,x)
  cf = return_cache(fa.f,x)
  cf, cg
end

function evaluate!(c,fa::CompositionFieldArrayField,x)#::AbstractArray{<:Point})
  cf, cg = c
  rg = evaluate!(cg,fa.g,x)
  rf = evaluate!(cf,fa.f,x)
end

# Gradients

# User API:
#
# BroadcastMapping(∇)(i_to_f)
#

@inline function evaluate!(cache,k::BroadcastMapping{typeof(gradient)},f::AbstractArray{<:Field})
  FieldGradientArray(f)
end

struct FieldGradientArray{A,T,N} <: AbstractArray{T,N}
  fa::A
  function FieldGradientArray(f::AbstractArray{<:Field})
    s = size(f)
    T = typeof(FieldGradient(eltype(f)))
    N = length(s)
    A = typeof(f)
    new{A,T,N}(f)
  end
end

@inline Base.size(a::FieldGradientArray) = size(a.fa)
@inline Base.axes(a::FieldGradientArray) = axes(a.fa)
@inline Base.getindex(a::FieldGradientArray,i::Integer...) = FieldGradient(a.fa[i...])
@inline Base.ndims(a::FieldGradientArray{A}) where A = ndims(A)
@inline Base.eltype(a::FieldGradientArray{A}) where A = FieldGradient{eltype(A)}
@inline Base.IndexStyle(::Type{<:FieldGradientArray{A}}) where A = IndexStyle(A)

@inline function evaluate!(cache,k::BroadcastMapping{typeof(gradient)},a::TransposeFieldVector)
  transpose(BroadcastMapping(∇)(a.basis))
end

@inline return_cache(f::FieldGradientArray,x) = return_gradient_cache(f.fa,x)

@inline evaluate!(cache,f::FieldGradientArray,x) = evaluate_gradient!(cache,f.fa,x)

@inline return_gradient_cache(fa::AbstractArray{<:Field},x) = return_cache(∇.(fa),x)

@inline evaluate_gradient!(cache,f::AbstractArray{<:Field},x) = evaluate!(cache,∇.(f),x)

# Hessian

struct FieldHessianArray{A,T,N} <: AbstractArray{T,N}
  fa::A
  function FieldHessianArray(f::AbstractArray{<:Field})
    s = size(f)
    T = typeof(FieldHessian(eltype(f)))
    N = length(s)
    A = typeof(f)
    new{A,T,N}(f)
  end
end

@inline return_hessian_cache(fa::AbstractArray{<:Field},x) = return_cache(∇.(∇.(fa)),x)

@inline evaluate_hessian!(cache,f::AbstractArray{<:Field},x) = evaluate!(cache,∇.(∇.(f)),x)


@inline Base.size(a::FieldHessianArray) = size(a.fa)
@inline Base.axes(a::FieldHessianArray) = axes(a.fa)
@inline Base.getindex(a::FieldHessianArray,i::Integer...) = FieldHessian(a.fa[i...])
@inline Base.ndims(a::FieldHessianArray{A}) where A = ndims(A)
@inline Base.eltype(a::FieldHessianArray{A}) where A = FieldHessian{eltype(A)}
@inline Base.IndexStyle(::Type{<:FieldHessianArray{A}}) where A = IndexStyle(A)

# @santiagobadia : Gradients of the previous operations ? needed?
# reimplement again chain rules etc for arrays... !

# Meeting 22 Sep
# integrate(f::Field,x::AbstractVector{<:Point},w::AbstractVector{<:Real}) = sum( f(x) .* w  )
# integrate(f::AbstractArray{<:Field},x::AbstractVector{<:Point},w::AbstractVector{<:Real}) =

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
