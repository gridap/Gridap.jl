
# Optimization for
#
#  g = lazy_map( linear_combination, cell_to_i_to_val,cell_to_i_to_f)
#  lazy_map(evaluate,g,cell_to_x)
#
function lazy_map(
  ::typeof(evaluate), a::LazyArray{<:Fill{typeof(linear_combination)}}, x::AbstractArray)

  i_to_values = a.f[1]
  i_to_basis = a.f[2]
  i_to_basis_x = lazy_map(evaluate,i_to_basis,x)
  lazy_map(LinearCombinationMap(:),i_to_values,i_to_basis_x)
end

#linear_combination_on_values(b_pi::AbstractMatrix,v_i::AbstractVector) = evaluate!(c,MatMul(),b_pi,v_i)
#linear_combination_on_values(b_pi::AbstractMatrix,v_ij::AbstractMatrix) = evaluate!(c,MatMul(),b_pi,v_ij)
#linear_combination_on_values(b_i::AbstractVector,v_i::AbstractVector) = b_i ⋅ v_i
#linear_combination_on_values(b_i::AbstractVector,v_ij::AbstractMatrix) = evaluate!(c,MatMul(),transpose(v_ij),b_i)


# Optimization for
#
#  g = lazy_map(transpose,cell_to_i_to_f)
#  lazy_map(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate), a::LazyArray{<:Fill{typeof(transpose)}}, x::AbstractArray)

  i_to_basis_x = lazy_map(evaluate,a.f[1],x)
  lazy_map(TransposeMap(),i_to_basis_x)
end

# Optimization for
#
#  g = lazy_map(∘,cell_to_f,cell_to_h)
#  lazy_map(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate), a::LazyArray{<:Fill{typeof(∘)}}, x::AbstractArray)

  f = a.f[1]
  g = a.f[2]
  gx = lazy_map(evaluate,g,x)
  fx = lazy_map(evaluate,f,gx)
  fx
end

# Optimization for
#
#  g = lazy_map(Broadcasting(∘),cell_to_i_to_f,cell_to_h)
#  lazy_map(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate), a::LazyArray{<:Fill{Broadcasting{typeof(∘)}}}, x::AbstractArray)

  f = a.f[1]
  g = a.f[2]
  gx = lazy_map(evaluate,g,x)
  fx = lazy_map(evaluate,f,gx)
  fx
end

# Optimization for
#
#  g = lazy_map(Operation(+),cell_to_f,cell_to_h)
#  lazy_map(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate),
  a::LazyArray{<:Fill{<:Operation}},
  x::AbstractArray)

  fx = map( fi->lazy_map(evaluate,fi,x), a.f)
  op = a.g.value.op
  lazy_map( Broadcasting(op), fx...)
end


# Optimization for
#
#  g = lazy_map(Broadcasting(Operation(+)),cell_to_f,cell_to_h)
#  lazy_map(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate), a::LazyArray{<:Fill{<:Broadcasting{<:Operation}}}, x::AbstractArray)

  fx = map( fi->lazy_map(evaluate,fi,x), a.f)
  op = a.g.value.f.op
  lazy_map(BroadcastingFieldOpMap(op),fx...)
end

# Optimization for
#
#  g = lazy_map( linear_combination, cell_to_i_to_val, cell_to_i_to_f)
#  lazy_map(gradient,g)
#
function lazy_map(
  ::typeof(gradient), a::LazyArray{<:Fill{typeof(linear_combination)}})

  i_to_basis = lazy_map(Broadcasting(∇),a.f[2])
  i_to_values = a.f[1]
  lazy_map(linear_combination,i_to_values,i_to_basis)
end

# Optimization for
#
#  g = lazy_map( linear_combination, cell_to_i_to_val, cell_to_i_to_f)
#  lazy_map(Broadcasting(gradient),g)
#
function lazy_map(
  k::Broadcasting{typeof(∇)}, a::LazyArray{<:Fill{typeof(linear_combination)}})

  i_to_basis = lazy_map(k,a.f[2])
  i_to_values = a.f[1]
  lazy_map(linear_combination,i_to_values,i_to_basis)
end

function lazy_map(
  k::Broadcasting{typeof(∇∇)}, a::LazyArray{<:Fill{typeof(linear_combination)}})

  i_to_basis = lazy_map(k,a.f[2])
  i_to_values = a.f[1]
  lazy_map(linear_combination,i_to_values,i_to_basis)
end

# Optimization for
#
#  g = lazy_map(transpose,cell_to_i_to_f)
#  lazy_map(Broadcasting(gradient),g)
#
function lazy_map(
  k::Broadcasting{typeof(∇)}, a::LazyArray{<:Fill{typeof(transpose)}})

  i_to_basis = lazy_map(k,a.f[1])
  lazy_map( transpose, i_to_basis)
end

function lazy_map(
  k::Broadcasting{typeof(∇∇)}, a::LazyArray{<:Fill{typeof(transpose)}})

  i_to_basis = lazy_map(k,a.f[1])
  lazy_map( transpose, i_to_basis)
end

# Gradient rules
for op in (:+,:-)
  @eval begin

    function lazy_map(
      ::typeof(gradient), a::LazyArray{<:Fill{Operation{typeof($op)}}})

      f = a.f
      g = map(i->lazy_map(gradient,i),f)
      lazy_map(Operation($op),g...)
    end

    function lazy_map(
      ::Broadcasting{typeof(gradient)}, a::LazyArray{<:Fill{Broadcasting{Operation{typeof($op)}}}})

      f = a.f
      g = map(i->lazy_map(Broadcasting(∇),i),f)
      lazy_map(Broadcasting(Operation($op)),g...)
    end

  end
end

for op in (:*,:⋅,:⊙,:⊗)
  @eval begin

    function lazy_map(
      ::typeof(gradient), a::LazyArray{<:Fill{Operation{typeof($op)}}})

      f = a.f
      @notimplementedif length(f) != 2
      g = map(i->lazy_map(gradient,i),f)
      k(F1,F2,G1,G2) = product_rule($op,F1,F2,G1,G2)
      lazy_map(Operation(k),f...,g...)
    end

    function lazy_map(
      ::Broadcasting{typeof(gradient)}, a::LazyArray{<:Fill{Broadcasting{Operation{typeof($op)}}}})

      f = a.f
      @notimplementedif length(f) != 2
      g = map(i->lazy_map(Broadcasting(∇),i),f)
      k(F1,F2,G1,G2) = product_rule($op,F1,F2,G1,G2)
      lazy_map(Broadcasting(Operation(k)),f...,g...)
    end

  end
end

function lazy_map(
  ::Broadcasting{typeof(gradient)}, a::LazyArray{<:Fill{Broadcasting{typeof(∘)}}})

  f = a.f[1]
  g = a.f[2]
  ∇f = lazy_map(Broadcasting(∇),f)
  h = lazy_map(Broadcasting(∘),∇f,g)
  ∇g = lazy_map(Broadcasting(∇),g)

  lazy_map(Broadcasting(Operation(⋅)),∇g,h)
end

function lazy_map(
  ::typeof(gradient), a::LazyArray{<:Fill{typeof(∘)}})

  f = a.f[1]
  g = a.f[2]
  ∇f = lazy_map(∇,f)
  h = lazy_map(∘,∇f,g)
  ∇g = lazy_map(∇,g)

  lazy_map(Operation(⋅),∇g,h)
end

# Integration

function lazy_map(
  ::typeof(integrate),f::AbstractArray,x::AbstractArray,w::AbstractArray)

  fx = lazy_map(evaluate,f,x)
  lazy_map(IntegrationMap(),fx,w)
end

function lazy_map(
  ::typeof(integrate),f::AbstractArray,x::AbstractArray,w::AbstractArray,j::AbstractArray)

  fx = lazy_map(evaluate,f,x)
  jx = lazy_map(evaluate,j,x)
  lazy_map(IntegrationMap(),fx,w,jx)
end

# Pushing the gradient

function lazy_map(k::typeof(push_∇),cell_∇a::AbstractArray,cell_map::AbstractArray)
  lazy_map(Broadcasting(push_∇),cell_∇a,cell_map)
end

function lazy_map(k::Broadcasting{typeof(push_∇)},cell_∇a::AbstractArray,cell_map::AbstractArray)
  cell_Jt = lazy_map(∇,cell_map)
  cell_invJt = lazy_map(Operation(inv),cell_Jt)
  lazy_map(Broadcasting(Operation(⋅)),cell_invJt,cell_∇a)
end

# Composing by the identity

function lazy_map(
  k::Broadcasting{typeof(∘)},
  ::Type{T},
  a::AbstractArray,
  b::Fill{<:GenericField{typeof(identity)}}) where T
  @assert length(a) == length(b)
  a
end

## Reindex
#
##  j_to_f = lazy_map(Reindex(i_to_f),j_to_i)
##  lazy_map(evaluate,j_to_f)
#function lazy_map(
#  ::typeof(evaluate), a::LazyArray{<:Fill{<:Reindex}},x::AbstractArray)
#
#  i_to_f = a.g.value.values
#  j_to_i = a.f[1]
#  i_to_fx = lazy_map(evaluate,i_to_f,i_to_x)
#  j_to_fx = lazy_map(Reindex(i_to_fx),j_to_i)
#  j_to_x = lazy_map(Reindex(i_to_fx),j_to_i)
#end

#function lazy_map(
#  ::Broadcasting{typeof(gradient)},
#  a::LazyArray{<:Fill{Broadcasting{Operation{typeof(*)}}}})
#
#  f = a.f
#  g = map(i->lazy_map(gradient,i),f)
#  r1 = lazy_map(Broadcasting(Operation(*)),f[1],g[2])
#  r2 = lazy_map(Broadcasting(Operation(*)),f[2],g[1])
#  lazy_map(Broadcasting(Operation(+)),r1,r2)
#end
#
#function lazy_map(
#  ::typeof(gradient), a::LazyArray{<:Fill{Operation{typeof(*)}}})
#
#  f = a.f
#  g = map(i->lazy_map(gradient,i),f)
#  r1 = lazy_map(Operation(*),f[1],g[2])
#  r2 = lazy_map(Operation(*),f[2],g[1])
#  lazy_map(Operation(+),r1,r2)
#end
#
## Function just used for dispatching
#integrate(f::Field,w,j,x) = transpose(evaluate(f,x))*(w.*meas.(j(x)))
#integrate(f::AbstractArray{<:Field},w,j,x) = transpose(evaluate(f,x))*(w.*meas.(j(x)))
#
## @santiagobadia : I would say that the integration points should be at our disposal
## when creating the LazyArray
#function lazy_map(
#  ::typeof(evaluate), a::LazyArray{<:Fill{typeof(integrate)}})#, x::AbstractArray)
#  f, w, j, x = a.f
#  fx = lazy_map(evaluate,f,x)
#  jx = lazy_map(evaluate,j,x)
#  k = Integrate()
#  lazy_map(k,fx,w,jx)
#end
#
## Other optimizations
##
#function lazy_map(
#  ::typeof(linear_combination), a::LazyArray{<:Fill{typeof(transpose)}}, i_to_values::AbstractArray)
#
#  i_to_basis = a.f[1]
#  lazy_map(linear_combination,i_to_basis,i_to_values)
#end
