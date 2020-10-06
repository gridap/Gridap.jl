
# optimization for
#
#    apply(evaluate,cell_to_f,cell_to_x)
#
# function apply(::typeof(evaluate),f::AbstractArray,x::AbstractArray)
  # apply(f,x)
# end

# Optimization for
#
#  g = apply( linear_combination, cell_to_i_to_f, cell_to_i_to_val)
#  apply(evaluate,g,cell_to_x)
#
function lazy_map(
  ::typeof(evaluate), a::MappedArray{<:Fill{typeof(linear_combination)}}, x::AbstractArray)

  i_to_basis = lazy_map(evaluate,a.f[1],x)
  i_to_values = a.f[2]
  lazy_map(LinCombVal(),i_to_basis,i_to_values)
end

linear_combination_on_values(b_pi::AbstractMatrix,v_i::AbstractVector) = evaluate!(c,MatMul(),b_pi,v_i)
linear_combination_on_values(b_pi::AbstractMatrix,v_ij::AbstractMatrix) = evaluate!(c,MatMul(),b_pi,v_ij)
linear_combination_on_values(b_i::AbstractVector,v_i::AbstractVector) = b_i ⋅ v_i
linear_combination_on_values(b_i::AbstractVector,v_ij::AbstractMatrix) = evaluate!(c,MatMul(),transpose(v_ij),b_i)


# Optimization for
#
#  g = apply(transpose,cell_to_i_to_f)
#  apply(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate), a::MappedArray{<:Fill{typeof(transpose)}}, x::AbstractArray)

  i_to_basis = lazy_map(evaluate,a.f[1],x)
  lazy_map(transpose_field_indices,i_to_basis)
end

# Optimization for
#
#  g = apply(∘,cell_to_f,cell_to_h)
#  apply(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate), a::MappedArray{<:Fill{typeof(∘)}}, x::AbstractArray)

  f = a.f[1]
  g = a.f[2]
  gx = lazy_map(evaluate,g,x)
  fx = lazy_map(evaluate,f,gx)
  fx
end

# Optimization for
#
#  g = apply(Broadcasting(∘),cell_to_i_to_f,cell_to_h)
#  apply(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate), a::MappedArray{<:Fill{Broadcasting{typeof(∘)}}}, x::AbstractArray)

  f = a.f[1]
  g = a.f[2]
  gx = lazy_map(evaluate,g,x)
  fx = lazy_map(evaluate,f,gx)
  fx
end

# Optimization for
#
#  g = apply(Operation(+),cell_to_f,cell_to_h)
#  apply(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate),
  a::MappedArray{<:Fill{<:Operation}},
  x::AbstractArray)

  fx = map( fi->lazy_map(evaluate,fi,x), a.f)
  op = a.g.value.op
  lazy_map( op, fx...)
end


# Optimization for
#
#  g = apply(Broadcasting(Operation(+)),cell_to_f,cell_to_h)
#  apply(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate), a::MappedArray{<:Fill{<:Broadcasting{<:Operation}}}, x::AbstractArray)

  fx = map( fi->lazy_map(evaluate,fi,x), a.f)
  op = Broadcasting(a.g.value.f.op)
  lazy_map(op,fx...)
end

# Optimization for
#
#  g = apply( linear_combination, cell_to_i_to_f, cell_to_i_to_val)
#  apply(gradient,g)
#
function lazy_map(
  ::typeof(gradient), a::MappedArray{<:Fill{typeof(linear_combination)}})

  i_to_basis = lazy_map(Broadcasting(gradient),a.f[1])
  i_to_values = a.f[2]
  lazy_map(linear_combination,i_to_basis,i_to_values)
end

# Optimization for
#
#  g = apply( linear_combination, cell_to_i_to_f, cell_to_i_to_val)
#  apply(Broadcasting(gradient),g)
#
function lazy_map(
  ::Broadcasting{typeof(gradient)}, a::MappedArray{<:Fill{typeof(linear_combination)}})

  i_to_basis = lazy_map(Broadcasting(gradient),a.f[1])
  i_to_values = a.f[2]
  lazy_map(linear_combination,i_to_basis,i_to_values)
end

# Optimization for
#
#  g = apply(transpose,cell_to_i_to_f)
#  apply(Broadcasting(gradient),g)
#
function lazy_map(
  ::Broadcasting{typeof(gradient)}, a::MappedArray{<:Fill{typeof(transpose)}})

  i_to_basis = lazy_map(gradient,a.f[1])
  lazy_map( transpose, i_to_basis)
end

# Product rules
for op in (:+,:-)
  @eval begin

    function lazy_map(
      ::typeof(gradient), a::MappedArray{<:Fill{Operation{typeof($op)}}})

      f = a.f
      g = map(i->lazy_map(gradient,i),f)
      lazy_map(Operation($op),g...)
    end

    function lazy_map(
      ::Broadcasting{typeof(gradient)}, a::MappedArray{<:Fill{Broadcasting{Operation{typeof($op)}}}})

      f = a.f
      g = map(i->lazy_map(gradient,i),f)
      lazy_map(Broadcasting(Operation($op)),g...)
    end

  end
end

function lazy_map(
  ::Broadcasting{typeof(gradient)},
  a::MappedArray{<:Fill{Broadcasting{Operation{typeof(*)}}}})

  f = a.f
  g = map(i->lazy_map(gradient,i),f)
  r1 = lazy_map(Broadcasting(Operation(*)),f[1],g[2])
  r2 = lazy_map(Broadcasting(Operation(*)),f[2],g[1])
  lazy_map(Broadcasting(Operation(+)),r1,r2)
end

function lazy_map(
  ::typeof(gradient), a::MappedArray{<:Fill{Operation{typeof(*)}}})

  f = a.f
  g = map(i->lazy_map(gradient,i),f)
  r1 = lazy_map(Operation(*),f[1],g[2])
  r2 = lazy_map(Operation(*),f[2],g[1])
  lazy_map(Operation(+),r1,r2)
end

# Function just used for dispatching
integrate(f::Field,w,j,x) = transpose(evaluate(f,x))*(w.*meas.(j(x)))
integrate(f::AbstractArray{<:Field},w,j,x) = transpose(evaluate(f,x))*(w.*meas.(j(x)))

# @santiagobadia : I would say that the integration points should be at our disposal
# when creating the MappedArray
function lazy_map(
  ::typeof(evaluate), a::MappedArray{<:Fill{typeof(integrate)}})#, x::AbstractArray)
  f, w, j, x = a.f
  fx = lazy_map(evaluate,f,x)
  jx = lazy_map(evaluate,j,x)
  k = Integrate()
  lazy_map(k,fx,w,jx)
end

# Other optimizations
#
function lazy_map(
  ::typeof(linear_combination), a::MappedArray{<:Fill{typeof(transpose)}}, i_to_values::AbstractArray)

  i_to_basis = a.f[1]
  lazy_map(linear_combination,i_to_basis,i_to_values)
end
