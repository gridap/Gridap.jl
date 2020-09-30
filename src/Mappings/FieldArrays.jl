# # Default implementation of evaluate! for arrays of fields

function return_cache(f::AbstractArray{<:Field},x::Point)
  T = return_type(first(f),first(x))
  c = CachedArray(zeros(T,size(f)))
end

function evaluate!(c,f::AbstractArray{<:Field},x::Point)
  setsize!(c,size(f))
  for j in eachindex(f)
    @inbounds c.array[j] = evaluate(f[j],x)
  end
  c.array
end

function return_cache(f::AbstractArray{<:Field},x::AbstractArray{<:Point})
  T = return_type(first(f),first(x))
  c = CachedArray(zeros(T,(size(x)...,size(f)...)))
end

function evaluate!(c,f::AbstractArray{<:Field},x::AbstractArray{<:Point})
  setsize!(c,(size(x)...,size(f)...))
  for j in eachindex(f)
    for i in eachindex(x)
      @inbounds c.array[i,j] = evaluate(f[j],x[i])
    end
  end
  c.array
end

# Broadcasting operations

# Non-performant default implementation

return_cache(b::BroadcastMapping,args::Union{Field,AbstractArray{<:Field}}...) = nothing

function evaluate!(cache,b::BroadcastMapping,args::Union{Field,AbstractArray{<:Field}}...)
  broadcast(b.f,args...)
end

# @santiagobadia : Here if we put AbstractArray it goes to the + in Julia for Array...
# for op in (:+,:-,:*,:⋅)
  # @eval $op(a::Array{<:Field},b::Array{<:Field}) = BroadcastMapping(Operation($op))(a,b)
# end

# Transpose

function Base.transpose(a::AbstractVector{<:Field})
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

Base.size(a::TransposeFieldVector) = (1,length(a.basis))
Base.axes(a::TransposeFieldVector) = (Base.OneTo(1),axes(a.basis,1))
Base.getindex(a::TransposeFieldVector,i::Integer,j::Integer) = a.basis[j,i]
Base.getindex(a::TransposeFieldVector,i::Integer) = a.basis[i]
Base.IndexStyle(::TransposeFieldVector{A}) where A = IndexStyle(A)

return_cache(a::TransposeFieldVector,x::Point) = return_cache(a.basis,x)
return_cache(a::TransposeFieldVector,x::AbstractArray{<:Point}) = return_cache(a.basis,x)

# @santiagobadia : How to fix this?
function evaluate!(cache,a::TransposeFieldVector,x::Point)
  M = evaluate!(cache,a.basis,x)
  transpose_field_indices(M)
end

function evaluate!(cache,a::TransposeFieldVector,x::AbstractArray{<:Point})
  M = evaluate!(cache,a.basis,x)
  transpose_field_indices(M)
end

transpose_field_indices(M::AbstractVector) = transpose(M)
transpose_field_indices(M::AbstractMatrix) = TransposeFieldIndices(M)

struct TransposeFieldIndices{A,T} <: AbstractArray{T,3}
  matrix::A
  @inline function TransposeFieldIndices(matrix::AbstractMatrix{T}) where T
    A = typeof(matrix)
    new{A,T}(matrix)
  end
end

Base.size(a::TransposeFieldIndices) = (size(a.matrix,1),1,size(a.matrix,2))
Base.axes(a::TransposeFieldIndices) = (axes(a.matrix,1),Base.OneTo(1),axes(a.matrix,2))
Base.IndexStyle(::Type{<:TransposeFieldIndices{A}}) where A = IndexStyle(A)
@inline Base.getindex(a::TransposeFieldIndices,i::Integer,j::Integer,k::Integer) = a.matrix[i,k]
@inline Base.getindex(a::TransposeFieldIndices,i::Integer) = a.matrix[i]
@inline Base.setindex!(a::TransposeFieldIndices,v,i::Integer,j::Integer,k::Integer) = (a.matrix[i,k] = v)
@inline Base.setindex!(a::TransposeFieldIndices,v,i::Integer) = (a.matrix[i] = v)

# Broadcast operations

(b::BroadcastMapping{<:Operation})(args::Union{Field,AbstractArray{<:Field}}...) = BroadcastOpFieldArray(b.f.op,args...)

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
# @santiagobadia : What do we want to do here?
@inline Base.IndexStyle(::Type{<:BroadcastOpFieldArray}) = IndexLinear
@inline Base.getindex(a::BroadcastOpFieldArray,I...) = broadcast(a.op,a.args...)[I...]

function return_cache(f::BroadcastOpFieldArray,x::Point)
  cfs = return_caches(f.args,x)
  rs = map((ci,fi) -> evaluate!(ci,fi,x),cfs,f.args)
  bm = BroadcastMapping(f.op)
  r = return_cache(bm,rs...)
  r, cfs
end

function evaluate!(c,f::BroadcastOpFieldArray,x::Point)
  r, cfs = c
  rs = map((ci,fi) -> evaluate!(ci,fi,x),cfs,f.args)
  bm = BroadcastMapping(f.op)
  evaluate!(r,bm,rs...)
  r.array
end

function return_cache(f::BroadcastOpFieldArray,x::AbstractArray{<:Point})
  cfs = return_caches(f.args,x)
  rs = map((ci,fi) -> evaluate!(ci,fi,x),cfs,f.args)
  bm = BroadcastMapping(f.op)
  r = return_cache(bm,rs...)
  r, cfs
end

function evaluate!(c,f::BroadcastOpFieldArray,x::AbstractArray{<:Point})
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

# @santiagobadia : How to implement nicely?
*(f::TransposeFieldVector,g::AbstractVector{<:Field}) = DotOpFieldVectors(f.basis,g)
⋅(f::TransposeFieldVector,g::AbstractVector{<:Field}) = DotOpFieldVectors(f.basis,g)
⋅(g::AbstractVector{<:Field},f::AbstractVector{<:Field}) = DotOpFieldVectors(f,g)
⋅(g::AbstractVector{<:Number},f::AbstractVector{<:Field}) = DotOpFieldVectors(f,GenericField.(g))
⋅(f::AbstractVector{<:Field},g::AbstractVector{<:Number}) = DotOpFieldVectors(f,GenericField.(g))

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

linear_combination(a,b) = a⋅b

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

# Do we need CompositionFieldFieldArray?
# function CompositionFieldArray(Field,AbstractArray{<:Field})

(b::BroadcastMapping{typeof(∘)})(f::AbstractArray{<:Field},g::Field) = CompositionFieldArrayField(f,g)


function return_cache(fa::CompositionFieldArrayField,x::Point)
  cg = return_cache(fa.g,x)
  rg = evaluate!(cg,fa.g,x)
  cf = return_cache(fa.f,x)
  cf, cg
end

function return_cache(fa::CompositionFieldArrayField,x::AbstractArray{<:Point})
  cg = return_cache(fa.g,x)
  rg = evaluate!(cg,fa.g,x)
  cf = return_cache(fa.f,x)
  cf, cg
end

function evaluate!(c,fa::CompositionFieldArrayField,x::Point)
  cf, cg = c
  rg = evaluate!(cg,fa.g,x)
  rf = evaluate!(cf,fa.f,x)
end

function evaluate!(c,fa::CompositionFieldArrayField,x::AbstractArray{<:Point})
  cf, cg = c
  rg = evaluate!(cg,fa.g,x)
  rf = evaluate!(cf,fa.f,x)
end

# Gradients

# User API:
#
# BroadcastMapping(∇)(i_to_f)
#

function evaluate!(cache,k::BroadcastMapping{typeof(gradient)},f::AbstractArray{<:Field})
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

Base.size(a::FieldGradientArray) = size(a.fa)
Base.axes(a::FieldGradientArray) = axes(a.fa)
Base.getindex(a::FieldGradientArray,i::Integer...) = FieldGradient(a.fa[i...])
Base.ndims(a::FieldGradientArray{A}) where A = ndims(A)
Base.eltype(a::FieldGradientArray{A}) where A = FieldGradient{eltype(A)}
Base.IndexStyle(::Type{<:FieldGradientArray{A}}) where A = IndexStyle(A)

function evaluate!(cache,k::BroadcastMapping{typeof(gradient)},a::TransposeFieldVector)
  transpose(BroadcastMapping(∇)(a.basis))
end

@inline return_cache(f::FieldGradientArray,x::Point) = return_gradient_cache(f.fa,x)
@inline return_cache(f::FieldGradientArray,x::AbstractArray{<:Point}) = return_gradient_cache(f.fa,x)

@inline evaluate!(cache,f::FieldGradientArray,x::Point) = evaluate_gradient!(cache,f.fa,x)
@inline evaluate!(cache,f::FieldGradientArray,x::AbstractArray{<:Point}) = evaluate_gradient!(cache,f.fa,x)

return_gradient_cache(fa::AbstractArray{<:Field},x::Point) = return_cache(∇.(fa),x)
return_hessian_cache(fa::AbstractArray{<:Field},x::Point) = return_cache(∇.(∇.(fa)),x)

@inline evaluate_gradient!(cache,f::AbstractArray{<:Field},x::Point) = evaluate(cache,∇.(f),x)
@inline evaluate_hessian!(cache,f::AbstractArray{<:Field},x::Point) = evaluate(cache,∇.(∇.(f)),x)

return_gradient_cache(fa::AbstractArray{<:Field},x::AbstractArray{<:Point}) = return_cache(∇.(fa),x)
return_hessian_cache(fa::AbstractArray{<:Field},x::AbstractArray{<:Point}) = return_cache(∇.(∇.(fa)),x)

@inline evaluate_gradient!(cache,f::AbstractArray{<:Field},x::AbstractArray{<:Point}) = evaluate(cache,∇.(f),x)
@inline evaluate_hessian!(cache,f::AbstractArray{<:Field},x::AbstractArray{<:Point}) = evaluate(cache,∇.(∇.(f)),x)

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

Base.size(a::FieldHessianArray) = size(a.fa)
Base.axes(a::FieldHessianArray) = axes(a.fa)
Base.getindex(a::FieldHessianArray,i::Integer...) = FieldHessian(a.fa[i...])
Base.ndims(a::FieldHessianArray{A}) where A = ndims(A)
Base.eltype(a::FieldHessianArray{A}) where A = FieldHessian{eltype(A)}
Base.IndexStyle(::Type{<:FieldHessianArray{A}}) where A = IndexStyle(A)

# @santiagobadia : Gradients of the previous operations ? needed?
# reimplement again chain rules etc for arrays... !
