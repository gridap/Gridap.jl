
# optimization for
#
#    apply(evaluate,cell_to_f,cell_to_x)
#
function apply(::typeof(evaluate),f::AbstractArray,x::AbstractArray)
  apply(f,x)
end

# Optimization for
#
#  g = apply( linear_combination, cell_to_i_to_f, cell_to_i_to_val)
#  apply(evaluate,g,cell_to_x)
#
function apply(
  ::typeof(evaluate), a::MappedArray{<:Fill{typeof(linear_combination)}}, x::AbstractArray)

  i_to_basis = apply(evaluate,a.f[1],x)
  i_to_values = a.f[2]
  apply(LinCombVal(),i_to_basis,i_to_values)
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
function apply(
  ::typeof(evaluate), a::MappedArray{<:Fill{typeof(transpose)}}, x::AbstractArray)

  i_to_basis = apply(evaluate,a.f[1],x)
  apply(transpose_field_indices,i_to_basis)
end

# Optimization for
#
#  g = apply(∘,cell_to_f,cell_to_h)
#  apply(evaluate,g)
#
function apply(
  ::typeof(evaluate), a::MappedArray{<:Fill{typeof(∘)}}, x::AbstractArray)

  f = a.f[1]
  g = a.f[2]
  gx = apply(evaluate,g,x)
  fx = apply(evaluate,f,gx)
  fx
end

# Optimization for
#
#  g = apply(BroadcastMapping(∘),cell_to_i_to_f,cell_to_h)
#  apply(evaluate,g)
#
function apply(
  ::typeof(evaluate), a::MappedArray{<:Fill{BroadcastMapping{typeof(∘)}}}, x::AbstractArray)

  f = a.f[1]
  g = a.f[2]
  gx = apply(evaluate,g,x)
  fx = apply(evaluate,f,gx)
  fx
end

# Optimization for
#
#  g = apply(Operation(+),cell_to_f,cell_to_h)
#  apply(evaluate,g)
#
function apply(
  ::typeof(evaluate),
  a::MappedArray{<:Fill{<:Operation}},
  x::AbstractArray)

  fx = map( fi->apply(evaluate,fi,x), a.f)
  op = a.g.value.op
  apply( op, fx...)
end


# Optimization for
#
#  g = apply(BroadcastMapping(Operation(+)),cell_to_f,cell_to_h)
#  apply(evaluate,g)
#
function apply(
  ::typeof(evaluate), a::MappedArray{<:Fill{<:BroadcastMapping{<:Operation}}}, x::AbstractArray)

  fx = map( fi->apply(evaluate,fi,x), a.f)
  op = BroadcastMapping(a.g.value.f.op)
  apply(op,fx...)
end

# Optimization for
#
#  g = apply( linear_combination, cell_to_i_to_f, cell_to_i_to_val)
#  apply(gradient,g)
#
function apply(
  ::typeof(gradient), a::MappedArray{<:Fill{typeof(linear_combination)}})

  i_to_basis = apply(BroadcastMapping(gradient),a.f[1])
  i_to_values = a.f[2]
  apply(linear_combination,i_to_basis,i_to_values)
end

# Optimization for
#
#  g = apply( linear_combination, cell_to_i_to_f, cell_to_i_to_val)
#  apply(BroadcastMapping(gradient),g)
#
function apply(
  ::BroadcastMapping{typeof(gradient)}, a::MappedArray{<:Fill{typeof(linear_combination)}})

  i_to_basis = apply(BroadcastMapping(gradient),a.f[1])
  i_to_values = a.f[2]
  apply(linear_combination,i_to_basis,i_to_values)
end

# Optimization for
#
#  g = apply(transpose,cell_to_i_to_f)
#  apply(BroadcastMapping(gradient),g)
#
function apply(
  ::BroadcastMapping{typeof(gradient)}, a::MappedArray{<:Fill{typeof(transpose)}})

  i_to_basis = apply(gradient,a.f[1])
  apply( transpose, i_to_basis)
end

# Product rules
for op in (:+,:-)
  @eval begin

    function apply(
      ::typeof(gradient), a::MappedArray{<:Fill{Operation{typeof($op)}}})

      f = a.f
      g = map(i->apply(gradient,i),f)
      apply(Operation($op),g...)
    end

    function apply(
      ::BroadcastMapping{typeof(gradient)}, a::MappedArray{<:Fill{BroadcastMapping{Operation{typeof($op)}}}})

      f = a.f
      g = map(i->apply(gradient,i),f)
      apply(BroadcastMapping(Operation($op)),g...)
    end

  end
end

function apply(
  ::BroadcastMapping{typeof(gradient)},
  a::MappedArray{<:Fill{BroadcastMapping{Operation{typeof(*)}}}})

  f = a.f
  g = map(i->apply(gradient,i),f)
  r1 = apply(BroadcastMapping(Operation(*)),f[1],g[2])
  r2 = apply(BroadcastMapping(Operation(*)),f[2],g[1])
  apply(BroadcastMapping(Operation(+)),r1,r2)
end

function apply(
  ::typeof(gradient), a::MappedArray{<:Fill{Operation{typeof(*)}}})

  f = a.f
  g = map(i->apply(gradient,i),f)
  r1 = apply(Operation(*),f[1],g[2])
  r2 = apply(Operation(*),f[2],g[1])
  apply(Operation(+),r1,r2)
end

# Function just used for dispatching
integrate(f::Field,w,j,x) = transpose(evaluate(f,x))*(w.*meas.(j(x)))
integrate(f::AbstractArray{<:Field},w,j,x) = transpose(evaluate(f,x))*(w.*meas.(j(x)))

# @santiagobadia : I would say that the integration points should be at our disposal
# when creating the MappedArray
function apply(
  ::typeof(evaluate), a::MappedArray{<:Fill{typeof(integrate)}})#, x::AbstractArray)
  f, w, j, x = a.f
  fx = apply(evaluate,f,x)
  jx = apply(evaluate,j,x)
  k = Integrate()
  apply(k,fx,w,jx)
end

# Other optimizations
#
function apply(
  ::typeof(linear_combination), a::MappedArray{<:Fill{typeof(transpose)}}, i_to_values::AbstractArray)

  i_to_basis = a.f[1]
  apply(linear_combination,i_to_basis,i_to_values)
end
