
# Make arrays of field behave like Maps

@inline function return_cache(f::AbstractArray{T},x::Point) where T<:Field
  S = return_type(first(f),x)
  cr = CachedArray(zeros(S,size(f)))
  if isconcretetype(T)
    cf = return_cache(first(f),x)
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
@inline function evaluate!(c,f::AbstractArray{T},x::Point) where T<:Field
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
  S = return_type(first(f),testitem(x))
  cr = CachedArray(zeros(S,(size(x)...,size(f)...)))
  if isconcretetype(T)
    cf = return_cache(first(f),testitem(x))
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
    for j in eachindex(f)
      for i in eachindex(x)
        @inbounds r[i,j] = evaluate!(cf,f[j],x[i])
      end
    end
  else
    for j in eachindex(f)
      for i in eachindex(x)
        @inbounds r[i,j] = evaluate(f[j],x[i])
      end
    end
  end
  r
end

function testargs(f::AbstractArray{T},x::Point) where T<:Field
  testargs(first(f),x)
end

function testargs(f::AbstractArray{T},x::AbstractArray{<:Point}) where T<:Field
  testargs(first(f),x)
end

# Default implementation of the hook methods for gradient and hessian

function return_gradient_cache(f::AbstractArray{T},x::Point) where T <:Field
  S = return_gradient_type(first(f),x)
  cr = CachedArray(zeros(S,size(f)))
  if isconcretetype(T)
    cf = return_gradient_cache(first(f),x)
  else
    cf = nothing
  end
  cr, cf
end

function evaluate_gradient!(cache,f::AbstractArray{T},x::Point) where T <:Field
  cr, cf = c
  setsize!(cr,size(f))
  r = cr.array
  if isconcretetype(T)
    for j in eachindex(f)
      @inbounds r[j] = evaluate_gradient!(cf,f[j],x)
    end
  else
    for j in eachindex(f)
      @inbounds r[j] = evaluate_gradient(f[j],x)
    end
  end
  r
end

function return_gradient_cache(f::AbstractArray{T},x::AbstractArray{<:Point}) where T <:Field
  S = return_gradient_type(first(f),testitem(x))
  cr = CachedArray(zeros(S,(size(x)...,size(f)...)))
  if isconcretetype(T)
    cf = return_gradient_cache(first(f),testitem(x))
  else
    cf = nothing
  end
  cr, cf
end

function evaluate_gradient!(cache,f::AbstractArray{T},x::AbstractArray{<:Point}) where T <:Field
  cr, cf = c
  setsize!(cr,(size(x)...,size(f)...))
  r = cr.array
  if isconcretetype(T)
    for j in eachindex(f)
      for i in eachindex(x)
        @inbounds r[i,j] = evaluate_gradient!(cf,f[j],x[i])
      end
    end
  else
    for j in eachindex(f)
      for i in eachindex(x)
        @inbounds r[i,j] = evaluate_gradient(f[j],x[i])
      end
    end
  end
  r
end

function return_hessian_cache(f::AbstractArray{T},x::Point) where T <:Field
end

function evaluate_hessian!(cache,f::AbstractArray{T},x::Point) where T <:Field
end

function return_hessian_cache(f::AbstractArray{T},x::AbstractArray{<:Point}) where T <:Field
end

function evaluate_hessian!(cache,f::AbstractArray{T},x::AbstractArray{<:Point}) where T <:Field
end




function test_field(f::AbstractArray{<:Field}, x, v, cmp=(==); grad=nothing, hess=nothing)
  test_mapping(f,(x,),v,cmp)
  if grad != nothing
    test_mapping(Broadcasting(gradient)(f),(x,),grad,cmp)
  end
  if hess != nothing
    test_mapping(Broadcasting(hessian)(f),(x,),hess,cmp)
  end
end

## Broadcasting operations
#
## Non-performant default implementation
#
#@inline return_cache(b::Broadcasting,args::Union{Field,AbstractArray{<:Field}}...) = nothing
#
#@inline function evaluate!(cache,b::Broadcasting,args::Union{Field,AbstractArray{<:Field}}...)
#  broadcast(b.f,args...)
#end
#
## @santiagobadia : I am not sure we want this syntatic sugar
## @santiagobadia : Here if we put AbstractArray it goes to the + in Julia for Array...
## for op in (:+,:-,:*,:⋅)
## @eval $op(a::Array{<:Field},b::Array{<:Field}) = Broadcasting(Operation($op))(a,b)
## @eval $op(a::AbstractArray{<:Field},b::AbstractArray{<:Field}) = Broadcasting(Operation($op))(a,b)
## end
#
## Transpose
#
#"""
#It returns the transpose of a field array, which is another field array.
#"""
#@inline function Base.transpose(a::AbstractVector{<:Field})
#  TransposeFieldVector(a)
#end
#
#struct TransposeFieldVector{B,T} <: AbstractArray{T,2}
#  basis::B
#  function TransposeFieldVector(basis)
#    @assert ndims(basis) == 1
#    T = eltype(basis)
#    new{typeof(basis),T}(basis)
#  end
#end
#
#@inline Base.size(a::TransposeFieldVector) = (1,length(a.basis))
#@inline Base.axes(a::TransposeFieldVector) = (Base.OneTo(1),axes(a.basis,1))
#@inline Base.getindex(a::TransposeFieldVector,i::Integer,j::Integer) = a.basis[j,i]
#@inline Base.getindex(a::TransposeFieldVector,i::Integer) = a.basis[i]
#@inline Base.IndexStyle(::TransposeFieldVector{A}) where A = IndexStyle(A)
#
#@inline return_cache(a::TransposeFieldVector,x) = return_cache(a.basis,x)
#
#function evaluate!(cache,a::TransposeFieldVector,x)
#  M = evaluate!(cache,a.basis,x)
#  transpose_field_indices(M)
#end
#
#@inline transpose_field_indices(M::AbstractVector) = transpose(M)
#@inline transpose_field_indices(M::AbstractMatrix) = TransposeFieldIndices(M)
#
## Broadcast operations
#
#@inline (b::Broadcasting{<:Operation})(args::Union{Field,AbstractArray{<:Field}}...) = BroadcastOpFieldArray(b.f.op,args...)
#
#"""
#Type that represents a broadcast operation over a set of `AbstractArray{<:Field}`.
#The result is a sub-type of `AbstractArray{<:Field}`
#"""
#struct BroadcastOpFieldArray{O,T,S<:Field,N} <: AbstractArray{S,N}
#  op::O
#  args::T
#  function BroadcastOpFieldArray(op,args::Union{Field,AbstractArray{<:Field}}...)
#    fs = map(first,args)
#    S = typeof(op(fs...))
#    s = map(size,args)
#    bs = Base.Broadcast.broadcast_shape(s...)
#    N = length(bs)
#    T = typeof(args)
#    O = typeof(op)
#    new{O,T,S,N}(op,args)
#  end
#end
#
#@inline Base.size(a::BroadcastOpFieldArray) = Base.Broadcast.broadcast_shape(map(size,a.args)...)
#@inline Base.IndexStyle(::Type{<:BroadcastOpFieldArray}) = IndexLinear
#@inline Base.getindex(a::BroadcastOpFieldArray,I...) = broadcast(a.op,a.args...)[I...]
#
#function return_cache(f::BroadcastOpFieldArray,x)
#  cfs = map(fi -> return_cache(fi,x),f.args)
#  rs = map((ci,fi) -> evaluate!(ci,fi,x),cfs,f.args)
#  bm = Broadcasting(f.op)
#  r = return_cache(bm,rs...)
#  r, cfs
#end
#
#function evaluate!(c,f::BroadcastOpFieldArray,x)
#  r, cfs = c
#  rs = map((ci,fi) -> evaluate!(ci,fi,x),cfs,f.args)
#  bm = Broadcasting(f.op)
#  evaluate!(r,bm,rs...)
#  r.array
#end
#
## Dot product vectors
#
#"""
#Type that represents the Julia dot-product for two vectors of fields.
#The result is a sub-type of `Field`
#"""
#struct DotOpFieldVectors{F,G} <: Field
#  f::F
#  g::G
#end
#
#@inline *(f::TransposeFieldVector,g::AbstractVector{<:Field}) = DotOpFieldVectors(f.basis,g)
#@inline ⋅(f::TransposeFieldVector,g::AbstractVector{<:Field}) = DotOpFieldVectors(f.basis,g)
#@inline ⋅(g::AbstractVector{<:Field},f::AbstractVector{<:Field}) = DotOpFieldVectors(f,g)
#@inline ⋅(g::AbstractVector{<:Number},f::AbstractVector{<:Field}) = DotOpFieldVectors(f,GenericField.(g))
#@inline ⋅(f::AbstractVector{<:Field},g::AbstractVector{<:Number}) = DotOpFieldVectors(f,GenericField.(g))
#
#function return_cache(f::DotOpFieldVectors,x::Point)
#  f1 = f.f
#  f2 = f.g
#  c1 = return_cache(f1,x)
#  c2 = return_cache(f2,x)
#  c1, c2
#end
#
#function evaluate!(c,f::DotOpFieldVectors,x::Point)
#  cfs = c
#  c1, c2 = cfs
#  f1 = f.f
#  f2 = f.g
#  r = cr.array
#  r1 = evaluate!(c1,f1,x)
#  r2 = evaluate!(c2,f2,x)
#  r1⋅r2
#end
#
#function return_cache(f::DotOpFieldVectors,x::AbstractArray{<:Point})
#  f1 = f.f
#  f2 = f.g
#  c1 = return_cache(f1,x)
#  c2 = return_cache(f2,x)
#  r1 = evaluate(f1,first(x))
#  r2 = evaluate(f2,first(x))
#  r = CachedArray(zeros(typeof(r1⋅r2),length(x)))
#  r, (c1, c2)
#end
#
#function evaluate!(c,f::DotOpFieldVectors,x::AbstractArray{<:Point})
#  cr, cfs = c
#  c1, c2 = cfs
#  f1 = f.f
#  f2 = f.g
#  setsize!(cr,size(x))
#  r = cr.array
#  r1 = evaluate!(c1,f1,x)
#  r2 = evaluate!(c2,f2,x)
#  for i in 1:length(x)
#    @inbounds r[i] = view(r1,i,:)⋅view(r2,i,:)
#  end
#  r
#end
#
## linear combination
#
## Function to be used for dispatch in LazyArray
#
#"""
#    linear_combination(a,b)
#It returns the linear combination of a vector of fields for a given vector of
#  coefficients. The result is a `LinearCombination` field.
#"""
#@inline linear_combination(a::AbstractVector{<:Field},b::AbstractArray{<:Number}) = transpose(a)*b
#@inline linear_combination(a::AbstractArray{<:Number},b::AbstractVector{<:Field}) = transpose(b)*a
#
## mat product
##
## The only user API needed
## transpose(i_to_f)*i_to_vals
## transpose(i_to_f)*ij_to_vals
##
#"""
#    a*b
#
#    Idem as `linear_combination(a,b)`
#"""
#@inline function Base.:*(a::TransposeFieldVector,b::AbstractVector{<:Number})
#  LinearCombinationField(a.basis,b)
#end
#
#@inline function Base.:*(a::TransposeFieldVector,b::AbstractMatrix{<:Number})
#  LinearCombinationField(a.basis,b)
#end
#
## These ones are not needed
#Base.:*(a::AbstractMatrix{<:Field},b::AbstractVector{<:Number}) = @notimplemented
#Base.:*(a::AbstractMatrix{<:Field},b::AbstractMatrix{<:Number}) = @notimplemented
#
#"""
#Sub-type of `Field` that represents the linear combination of a vector of fields
#for a given vector of coefficients.
#"""
#struct LinearCombinationField{B,V} <: Field
#  basis::B
#  values::V
#end
#
#function return_cache(a::LinearCombinationField,x)
#  cb = return_cache(a.basis,x)
#  bx = evaluate!(cb, a.basis, x)
#  v = a.values
#  cr = return_cache(LinCombVal(),bx,v)
#  cb, cr
#end
#
#function evaluate!(cache,a::LinearCombinationField,x)
#  cb, cr = cache
#  bx = evaluate!(cb,a.basis,x)
#  evaluate!(cr,LinCombVal(),bx,a.values)
#end
#
## Composition
#
#"""
#Sub-type of `Field` that represents the composition of two fields. It can
#also represent the vector of fields that results from the composition of
#every field in a vector of fields and another field.
#"""
#struct CompositionFieldArrayField{T,N,A,B} <:AbstractArray{T,N}
#  f::A
#  g::B
#  function CompositionFieldArrayField(f::AbstractArray{<:Field},g::Field)
#    s = size(f)
#    T = typeof(first(f)∘g)
#    N = length(s)
#    A = typeof(f)
#    B = typeof(g)
#    new{T,N,A,B}(f,g)
#  end
#end
#
#"""
#Composition of a field (or vector of fields) and another field. It returns a
#`CompositionFieldArrayField`
#"""
#@inline (b::Broadcasting{typeof(∘)})(f::AbstractArray{<:Field},g::Field) = CompositionFieldArrayField(f,g)
#
#function return_cache(fa::CompositionFieldArrayField,x)
#  cg = return_cache(fa.g,x)
#  rg = evaluate!(cg,fa.g,x)
#  cf = return_cache(fa.f,x)
#  cf, cg
#end
#
#function evaluate!(c,fa::CompositionFieldArrayField,x)
#  cf, cg = c
#  rg = evaluate!(cg,fa.g,x)
#  rf = evaluate!(cf,fa.f,x)
#end
#
## Gradients
#
## User API:
##
## Broadcasting(∇)(i_to_f)
##
#
#@inline function evaluate!(cache,k::Broadcasting{typeof(gradient)},f::AbstractArray{<:Field})
#  FieldGradientArray(f)
#end
#
#"""
#A wrapper that represents the broadcast of `gradient` over an array of fields.
#"""
#struct FieldGradientArray{A,T,N} <: AbstractArray{T,N}
#  fa::A
#  function FieldGradientArray(f::AbstractArray{<:Field})
#    s = size(f)
#    T = typeof(FieldGradient(eltype(f)))
#    N = length(s)
#    A = typeof(f)
#    new{A,T,N}(f)
#  end
#end
#
#@inline Base.size(a::FieldGradientArray) = size(a.fa)
#@inline Base.axes(a::FieldGradientArray) = axes(a.fa)
#@inline Base.getindex(a::FieldGradientArray,i::Integer...) = FieldGradient(a.fa[i...])
#@inline Base.ndims(a::FieldGradientArray{A}) where A = ndims(A)
#@inline Base.eltype(a::FieldGradientArray{A}) where A = FieldGradient{eltype(A)}
#@inline Base.IndexStyle(::Type{<:FieldGradientArray{A}}) where A = IndexStyle(A)
#
#@inline function evaluate!(cache,k::Broadcasting{typeof(gradient)},a::TransposeFieldVector)
#  transpose(Broadcasting(∇)(a.basis))
#end
#
#@inline gradient(f::AbstractArray{<:Field}) = FieldGradientArray(f)
#
#@inline return_cache(f::FieldGradientArray,x) = return_gradient_cache(f.fa,x)
#
#@inline evaluate!(cache,f::FieldGradientArray,x) = evaluate_gradient!(cache,f.fa,x)
#
#@inline return_gradient_cache(fa::AbstractArray{<:Field},x) = return_cache(∇.(fa),x)
#
#@inline evaluate_gradient!(cache,f::AbstractArray{<:Field},x) = evaluate!(cache,∇.(f),x)
#
## Hessian
#
#"""
#A wrapper that represents the application of `gradient` twice
#over an array of fields.
#"""
#struct FieldHessianArray{A,T,N} <: AbstractArray{T,N}
#  fa::A
#  function FieldHessianArray(f::AbstractArray{<:Field})
#    s = size(f)
#    T = typeof(FieldHessian(eltype(f)))
#    N = length(s)
#    A = typeof(f)
#    new{A,T,N}(f)
#  end
#end
#
#@inline gradient(f::FieldGradientArray) = FieldHessianArray(f.fa)
#
#@inline hessian(f::AbstractArray{<:Field}) = FieldHessianArray(f)
#
#@inline return_cache(f::FieldHessianArray,x) = return_hessian_cache(f.fa,x)
#
#@inline evaluate!(cache,f::FieldHessianArray,x) = evaluate_hessian!(cache,f.fa,x)
#
#@inline return_hessian_cache(fa::AbstractArray{<:Field},x) = return_cache(∇.(∇.(fa)),x)
#
#@inline evaluate_hessian!(cache,f::AbstractArray{<:Field},x) = evaluate!(cache,∇.(∇.(f)),x)
#
#
#@inline Base.size(a::FieldHessianArray) = size(a.fa)
#@inline Base.axes(a::FieldHessianArray) = axes(a.fa)
#@inline Base.getindex(a::FieldHessianArray,i::Integer...) = FieldHessian(a.fa[i...])
#@inline Base.ndims(a::FieldHessianArray{A}) where A = ndims(A)
#@inline Base.eltype(a::FieldHessianArray{A}) where A = FieldHessian{eltype(A)}
#@inline Base.IndexStyle(::Type{<:FieldHessianArray{A}}) where A = IndexStyle(A)
#
## @santiagobadia : Gradients of the previous operations ? needed? we can leave it
## for the future if really needed
#
## Testers
#
#function test_field_array(f,p,nf;grad=false,hessian=false)
#  fa = fill(f,nf...)
#  fp = evaluate(f,p)
#  fap = fill(fp,nf...)
#  test_field(fa,p,fap)
#  if hessian
#    test_field_array(∇(f),p,nf;grad=true)
#  elseif grad
#    test_field_array(∇(f),p,nf)
#  end
#  fa, p
#end
#
#function test_broadcast_field_array(f,p,nf,np;grad=false,hessian=false)
#  test_field_array(f,p,nf,grad=grad,hessian=hessian)
#  fa = fill(f,nf...)
#  fp = evaluate(f,p)
#  x = fill(p,np...)
#  fax = fill(fp,np...,nf...)
#  test_field(fa,x,fax)
#  if hessian
#    test_broadcast_field_array(∇(f),p,nf,np;grad=true)
#  elseif grad
#    test_broadcast_field_array(∇(f),p,nf,np)
#  end
#  fa, x
#end
#
#function test_operation_field_array(op,x,fs...)
#  fop = Operation(op)
#  fb = op(fs...)
#  fba = fop(fs...)
#  if fb isa AbstractArray
#    @test eltype(fb) <: Field
#    @test fb isa AbstractArray{<:Field}
#    @test fba isa OperationArray
#    @test fba.res == fb
#  else
#    @test fb isa Field
#  end
#  @test evaluate(fb,x) == evaluate(fba,x)
#  fba
#end
#
#"""
#Given a matrix `np` x `nf1` x `nf2` result of the evaluation of a field vector
#on a vector of points, it returns an array in which the field axes (second and
#third axes) are permuted. It is equivalent as `Base.permutedims(A,(1,3,2)`
#but more performant, since it does not involve allocations.
#"""
#struct TransposeFieldIndices{A,T} <: AbstractArray{T,3}
#  matrix::A
#  @inline function TransposeFieldIndices(matrix::AbstractMatrix{T}) where T
#    A = typeof(matrix)
#    new{A,T}(matrix)
#  end
#end
#
#@inline Base.size(a::TransposeFieldIndices) = (size(a.matrix,1),1,size(a.matrix,2))
#@inline Base.axes(a::TransposeFieldIndices) = (axes(a.matrix,1),Base.OneTo(1),axes(a.matrix,2))
#@inline Base.IndexStyle(::Type{<:TransposeFieldIndices{A}}) where A = IndexStyle(A)
#@inline Base.getindex(a::TransposeFieldIndices,i::Integer,j::Integer,k::Integer) = a.matrix[i,k]
#@inline Base.getindex(a::TransposeFieldIndices,i::Integer) = a.matrix[i]
#@inline Base.setindex!(a::TransposeFieldIndices,v,i::Integer,j::Integer,k::Integer) = (a.matrix[i,k] = v)
#@inline Base.setindex!(a::TransposeFieldIndices,v,i::Integer) = (a.matrix[i] = v)
