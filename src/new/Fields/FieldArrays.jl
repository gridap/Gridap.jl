
"""
    evaluate_field_array(f::AbstractArray,x::AbstractArray) -> AbstractArray

Evaluates the fields in the array `f` at all the vector of points in the 
array of vector of points `x` and returns the result as a lazy array.

The result is numerically equivalent to 

    map(evaluate_field,a,x)
"""
function evaluate_field_array(a::AbstractArray,x::AbstractArray)
  _evaluate_field_array(a,x)
end

function _evaluate_field_array(a::AbstractArray,x::AbstractArray)
  k = Eval()
  apply(k,a,x)
end

struct Eval <: Kernel end

function kernel_cache(k::Eval,a,x)
  field_cache(a,x)
end

function apply_kernel!(cache,k::Eval,a,x)
  evaluate_field!(cache,a,x)
end

function kernel_return_type(k::Eval,a,x)
  field_return_type(a,x)
end

function _evaluate_field_array(a::AbstractArray{<:Field},x::AbstractArray)
  apply(a,x)
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
  a::AppliedArray{T,N,F,<:Fill} where {T,N,F},x::AbstractArray)
  kernel_evaluate(a.g.value,x,a.f...)
end

"""
    kernel_evaluate(k,x,f...)

Function to control the evaluation of the field resulting from the operation `apply(k,f...)`, where `k`
is a kernel and `f...` contains several arrays of fields.

By default, this function is implemented as

    g = apply(k,f...)
    evaluate_field_array(g,x)

However, it can be rewritten for specific kernels in order to improve performance and simplify
the underlying operation tree.

"""
function kernel_evaluate(k,x,f...)
  a = apply(k,f...)
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
  apply(k,a)
end


struct Grad <: Kernel end

@inline apply_kernel!(::Nothing,k::Grad,x) = field_gradient(x)

function field_array_gradient(
  a::AppliedArray{T,N,F,<:Fill} where {T,N,F})
  apply_gradient(a.g.value,a.f...)
end

"""
    apply_gradient(k,f...)

Function to control the generation of the gradient the field resulting from the operation `apply(k,f...)`, where `k`
is a kernel and `f...` contains several arrays of fields.

By default, it returns the array obtained as

    a = apply(k,f...)
    field_array_gradient(a)

However, it can be rewritten for specific kernels in order to improve performance and simplify
the underlying operation tree.

"""
function apply_gradient(k,f...)
  a = apply(k,f...)
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
  cfi = field_cache(fi,xi)
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
    fxi = evaluate_field!(cfi,fi,xi)
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
    apply_to_field_array(k,f::AbstractArray...)
    apply_to_field_array(::Type{T},k,f::AbstractArray...) where T

Returns an array of fields numerically equivalent to

    map( (x...) -> apply_kernel_to_field(k,x...), f )

"""
function apply_to_field_array(
  k,f::AbstractArray...)
  v = Valued(k)
  apply(v,f...)
end

function apply_to_field_array(
  ::Type{T},k,f::AbstractArray...) where T
  v = Valued(k)
  apply(T,v,f...)
end

struct Valued{K} <: Kernel
  k::K
  function Valued(k)
    new{typeof(k)}(k)
  end
end

@inline function apply_kernel!(cache,k::Valued,x...)
  apply_kernel_to_field(k.k,x...)
end

function kernel_evaluate(k::Valued,x,f...)
  fx = evaluate_field_arrays(f,x)
  a = apply(k.k,fx...)
end

for op in (:+,:-)
  @eval begin
    function apply_gradient(k::Valued{BCasted{typeof($op)}},f...)
      g = field_array_gradients(f...)
      apply(k,g...)
    end
  end
end

