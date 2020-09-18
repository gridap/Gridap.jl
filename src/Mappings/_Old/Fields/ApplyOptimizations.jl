# efficiently evaluate Fields in lazy arrays
# apply(evaluate,f,x)
# is chosen as the API that is optimized

function apply(::typeof(evaluate),f::AbstractArray{<:Field},x::AbstractArray)
  apply(f,x)
end


function apply(
  ::typeof(evaluate),
  a::MappedArray{<:Fill{typeof(linear_combination)}},
  x::AbstractArray)

  i_to_basis = apply(evaluate,a.f[1],x)
  i_to_values = a.f[2]
  apply(linear_combination_on_values,i_to_basis,i_to_values)
end

function apply(
  ::typeof(evaluate), a::MappedArray{<:Fill{typeof(trialize)}}, x::AbstractArray)

  i_to_basis = apply(evaluate,a.f[1],x)
  apply(trialize_on_values,i_to_basis)
end

function apply(
  ::typeof(evaluate),
  a::MappedArray{<:Fill{typeof(∘)}},
  x::AbstractArray)

  f = a.f[1]
  g = a.f[2]
  gx = apply(evaluate,g,x)
  fx = apply(evaluate,f,gx)
  fx
end

function apply(
  ::typeof(evaluate),
  a::MappedArray{<:Fill{<:FieldOperation}},
  x::AbstractArray)

  fx = map( fi->apply(evaluate,fi,x), a.f)
  op = a.g.value.op
  apply( (args...) -> broadcast(op,args...), fx... )  # TODO use BroadcastMapping
end

# Efficient gradients for Fields in lazy arrays
# apply(gradient,f)
# is the optimized API

function apply(
  ::typeof(gradient),
  a::MappedArray{<:Fill{typeof(linear_combination)}})
  i_to_basis = apply(gradient,a.f[1])
  i_to_values = a.f[2]
  apply(linear_combination,i_to_basis,i_to_values)
end

# TODO product rules
function apply(
  ::typeof(gradient), a::MappedArray{<:Fill{typeof(trialize)}})

  i_to_basis = apply(gradient,a.f[1])
  apply(trialize,i_to_basis)
end

# Product rules

function apply(
  ::typeof(gradient), a::MappedArray{<:Fill{FieldOperation(typeof(+))}})

  f = a.f
  g = map(i->apply(gradient,i),f)
  apply(FieldOperation(+),g...)
end

function apply(
  ::typeof(gradient), a::MappedArray{<:Fill{FieldOperation(typeof(-))}})

  f = a.f
  g = map(i->apply(gradient,i),f)
  apply(FieldOperation(-),g...)
end

function apply(
  ::typeof(gradient), a::MappedArray{<:Fill{FieldOperation(typeof(*))}})

  f = a.f
  g = map(i->apply(gradient,i),f)
  r1 = apply(FieldOperation(*),f[1],g[2])
  r2 = apply(FieldOperation(*),f[2],g[1])
  apply(FieldOperation(*),r1,r2)
end

# Other optimizations

function apply(
  ::typeof(linear_combination), a::MappedArray{<:Fill{typeof(trialize)}}, i_to_values::AbstractArray)

  i_to_basis = a.f[1]
  apply(linear_combination,i_to_basis,i_to_values)
end
