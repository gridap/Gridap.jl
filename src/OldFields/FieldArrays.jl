"""
    evaluate_field_array(f::AbstractArray,x::AbstractArray) -> AbstractArray

Evaluates the fields in the array `f` at all the vector of points in the
array of vector of points `x` and returns the result as a lazy array.

The result is numerically equivalent to

    map(evaluate_field,f,x)
"""
function evaluate_field_array(a::AbstractArray,x::AbstractArray)
  _evaluate_field_array(a,x)
end

function _evaluate_field_array(a::AbstractArray,x::AbstractArray)
  k = Eval()
  lazy_map(k,a,x)
end

struct Eval <: Mapping end

function return_cache(k::Eval,a,x)
  return_cache(a,x)
end

function evaluate!(cache,k::Eval,a,x)
  evaluate!(cache,a,x)
end

function return_type(k::Eval,a,x)
  field_return_type(a,x)
end

function _evaluate_field_array(a::AbstractArray{<:Field},x::AbstractArray)
  lazy_map(a,x)
end

"""
    evaluate_field_arrays(f::Tuple,x::AbstractArray) -> Tuple

Equivalent to

    tuple((evaluate_field_array(fi,x) for fi in f)...)
"""
function evaluate_field_arrays(f::Tuple,x::AbstractArray)
  _evaluate_field_arrays(x,f...)
end

function _evaluate_field_arrays(x,a,b...)
  ax = evaluate_field_array(a,x)
  bx = evaluate_field_arrays(b,x)
  (ax,bx...)
end

function _evaluate_field_arrays(x,a)
  ax = evaluate_field_array(a,x)
  (ax,)
end

"""
    evaluate(a::AbstractArray{<:Field},x::AbstractArray)

Equivalent to

    evaluate_field_array(a,x)

But only for arrays `a` whose element type inherits from `Field`. If this is not the case,
use `evaluate_field_array(a,x)` instead.
"""
function evaluate(a::AbstractArray{<:Field},x::AbstractArray)
  evaluate_field_array(a,x)
end

# Optimized version for arrays of fields obtained from a kernel
# and other arrays
function evaluate_field_array(
  a::LazyArray{T,N,F,<:Fill} where {T,N,F},x::AbstractArray)
  kernel_evaluate(a.g.value,x,a.f...)
end

"""
    kernel_evaluate(k,x,f...)

Function to control the evaluation of the field resulting from the operation `lazy_map(k,f...)`, where `k`
is a kernel and `f...` contains several arrays of fields.

By default, this function is implemented as

    g = lazy_map(k,f...)
    evaluate_field_array(g,x)

However, it can be rewritten for specific kernels in order to improve performance and simplify
the underlying operation tree.

"""
function kernel_evaluate(k,x,f...)
  a = lazy_map(k,f...)
  _evaluate_field_array(a,x)
end

"""
    field_array_gradient(a::AbstractArray)

Returns an array containing the gradients of the fields in the array `a`.
Numerically equivalent to

    map(field_gradient,a)
"""
function field_array_gradient(a::AbstractArray)
  _field_array_gradient(a)
end

function _field_array_gradient(a::AbstractArray)
  k = Grad()
  lazy_map(k,a)
end


struct Grad <: Mapping end

@inline evaluate!(::Nothing,k::Grad,x) = field_gradient(x)

function field_array_gradient(
  a::LazyArray{T,N,F,<:Fill} where {T,N,F})
  lazy_map_gradient(a.g.value,a.f...)
end

"""
    lazy_map_gradient(k,f...)

Function to control the generation of the gradient the field resulting from the operation `lazy_map(k,f...)`, where `k`
is a kernel and `f...` contains several arrays of fields.

By default, it returns the array obtained as

    a = lazy_map(k,f...)
    field_array_gradient(a)

However, it can be rewritten for specific kernels in order to improve performance and simplify
the underlying operation tree.

"""
function lazy_map_gradient(k,f...)
  a = lazy_map(k,f...)
  _field_array_gradient(a)
end

"""
    gradient(f::AbstractArray{<:Field})

Equivalent to

    field_array_gradient(f)

but only for arrays whose element type is `<:Field`. Use function `field_array_gradient` otherwise.
"""
function gradient(f::AbstractArray{<:Field})
  field_array_gradient(f)
end

"""
    field_array_gradients(f...)

Equivalent to

    map(field_array_gradient,f)
"""
function field_array_gradients(f...)
  map(field_array_gradient,f)
end

"""
    field_array_cache(a::AbstractArray,x::AbstractArray) -> Tuple

Returns the caches needed to perform the following iteration

    ca, cfi, cx = field_array_cache(a,x)

    for i in length(a)
      fi = getindex!(ca,a,i)
      xi = getindex!(cx,x,i)
      fxi = evaluate!(cfi,fi,xi)
    end
"""
function field_array_cache(a::AbstractArray,x::AbstractArray)
  ca = array_cache(a)
  fi = testitem(a)
  xi = testitem(x)
  cfi = return_cache(fi,xi)
  cx = array_cache(x)
  (ca,cfi,cx)
end

"""
    function test_array_of_fields(
      a::AbstractArray,
      x::AbstractArray,
      v::AbstractArray,
      cmp::Function=(==);
      grad = nothing)

Function to test an array of fields `a`. The array `v` is the expected result when calling
`evaluate_field_array(a,x)`. The entries in the computed array and the expected one are compared
with the `cmp` function. The key-word argument `grad` is optional. If present, it should contain
the expected result of

    ∇a = field_array_gradient(a)
    evaluate_field_array(∇a,x)
"""
function test_array_of_fields(
  a::AbstractArray,
  x::AbstractArray,
  v::AbstractArray,
  cmp::Function=(==);
  grad = nothing)

  ax = evaluate_field_array(a,x)
  test_array(ax,v,cmp)

  ca, cfi, cx = field_array_cache(a,x)

  t = true
  for i in 1:length(a)
    fi = getindex!(ca,a,i)
    xi = getindex!(cx,x,i)
    fxi = evaluate!(cfi,fi,xi)
    vi = v[i]
    ti = cmp(fxi,vi)
    t = t && ti
  end
  @test t

  if grad != nothing
    g = field_array_gradient(a)
    test_array_of_fields(g,x,grad,cmp)
  end

end

"""
    lazy_map_to_field_array(k,f::AbstractArray...)
    lazy_map_to_field_array(::Type{T},k,f::AbstractArray...) where T

Returns an array of fields numerically equivalent to

    map( (x...) -> evaluate_to_field(k,x...), f )

"""
function lazy_map_to_field_array(
  k,f::AbstractArray...)
  v = Valued(k)
  lazy_map(v,f...)
end

function lazy_map_to_field_array(
  ::Type{T},k,f::AbstractArray...) where T
  v = Valued(k)
  lazy_map(T,v,f...)
end

struct Valued{K} <: Mapping
  k::K
  function Valued(k)
    new{typeof(k)}(k)
  end
end

@inline function evaluate!(cache,k::Valued,x...)
  evaluate_to_field(k.k,x...)
end

function kernel_evaluate(k::Valued,x,f...)
  fx = evaluate_field_arrays(f,x)
  a = lazy_map(k.k,fx...)
end

# More optimizations

function evaluate_field_array(a::AppendedArray,b::AppendedArray)
  if (length(a.a) == length(b.a)) && (length(a.b) == length(b.b))
    c_a = evaluate_field_array(a.a,b.a)
    c_b = evaluate_field_array(a.b,b.b)
    lazy_append(c_a,c_b)
  else
    _evaluate_field_array(a,b)
  end
end

function evaluate_field_array(a::AppendedArray,b::AbstractArray)
  n = length(a.a)
  _b = lazy_append(lazy_split(b,n)...)
  evaluate_field_array(a,_b)
end

function field_array_gradient(a::AppendedArray)
  c_a = field_array_gradient(a.a)
  c_b = field_array_gradient(a.b)
  lazy_append(c_a,c_b)
end
