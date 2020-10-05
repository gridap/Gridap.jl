# Extended Array interface

function return_cache(a::AbstractArray,i...)
  nothing
end

getindex!(cache,a::AbstractArray,i...) = a[i...]

function testitem(a::AbstractArray{T}) where T
  if length(a) >0
    first(a)
  else
    testvalue(T)
  end::T
end

get_array(a::AbstractArray) = a

# Working with several arrays

@inline function getitems!(cf::Tuple,a::Tuple,i...)
  _getitems!(cf,i,a...)
end

getitems!(::Tuple{},::Tuple{},i) = ()

@inline function _getitems!(c,i,a,b...)
  ca = c[1]
  cb = c[2,end]
  ai = getindex!(ca,a,i...)
  bi = getitems!(cb,b,i...)
  (ai,bi...)
end

@inline function _getitems!(c,i,a)
  ca, = c
  ai = getindex!(ca,a,i...)
  (ai,)
end

@inline function getitems(a::Tuple,i...)
  _getitems(i,a...)
end

@inline function _getitems(i,a,b...)
  ai = a[i...]
  bi = getitems(b,i...)
  (ai,bi...)
end

@inline function _getitems(i,a)
  ai = a[i...]
  (ai,)
end
