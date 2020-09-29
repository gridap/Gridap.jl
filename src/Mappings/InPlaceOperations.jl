
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
  for i in 1:size(a)[1]
    _a = selectdim(a,1,i)
    _b = selectdim(a,1,i)
    _r = selectdim(r,1,i)
    mul!(_r,_a,_b)
  end
end

@inline function _inplace!(::typeof(*),r::AbstractArray,a::AbstractArray,b::AbstractArray,c::AbstractArray...)
  @notimplemented
end
