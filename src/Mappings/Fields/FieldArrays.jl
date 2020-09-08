


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

  ca, cfi, cx = return_mapping_array_cache(a,x)

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
