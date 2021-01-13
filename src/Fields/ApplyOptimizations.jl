
# Optimization for
#
#  g = lazy_map( linear_combination, cell_to_i_to_val,cell_to_i_to_f)
#  lazy_map(evaluate,g,cell_to_x)
#
function lazy_map(
  ::typeof(evaluate), a::LazyArray{<:Fill{typeof(linear_combination)}}, x::AbstractArray)

  i_to_values = a.args[1]
  i_to_basis = a.args[2]
  i_to_basis_x = lazy_map(evaluate,i_to_basis,x)
  lazy_map(LinearCombinationMap(:),i_to_values,i_to_basis_x)
end

# We always keep the parent when transposing (needed for some optimizations below)

function lazy_map(k::typeof(transpose),f::AbstractArray)
  fi = testitem(f)
  T = return_type(k,fi)
  s = size(f)
  LazyArray(T,Fill(k,s),f)
end

# Optimization for
#
#  g = lazy_map(transpose,cell_to_i_to_f)
#  lazy_map(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate), a::LazyArray{<:Fill{typeof(transpose)}}, x::AbstractArray)

  i_to_basis_x = lazy_map(evaluate,a.args[1],x)
  lazy_map(TransposeMap(),i_to_basis_x)
end

# Optimization for
#
#  g = lazy_map(∘,cell_to_f,cell_to_h)
#  lazy_map(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate), a::LazyArray{<:Fill{typeof(∘)}}, x::AbstractArray)

  f = a.args[1]
  g = a.args[2]
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

  f = a.args[1]
  g = a.args[2]
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

  fx = map( fi->lazy_map(evaluate,fi,x), a.args)
  op = a.maps.value.op
  lazy_map( Broadcasting(op), fx...)
end

# Optimization for
#
#  g = lazy_map(Broadcasting(Operation(+)),cell_to_f,cell_to_h)
#  lazy_map(evaluate,g)
#
function lazy_map(
  ::typeof(evaluate), a::LazyArray{<:Fill{<:Broadcasting{<:Operation}}}, x::AbstractArray)

  fx = map( fi->lazy_map(evaluate,fi,x), a.args)
  op = a.maps.value.f.op
  lazy_map(BroadcastingFieldOpMap(op),fx...)
end

# Optimization for
#
#  g = lazy_map( linear_combination, cell_to_i_to_val, cell_to_i_to_f)
#  lazy_map(gradient,g)
#
function lazy_map(
  ::typeof(gradient), a::LazyArray{<:Fill{typeof(linear_combination)}})

  i_to_basis = lazy_map(Broadcasting(∇),a.args[2])
  i_to_values = a.args[1]
  lazy_map(linear_combination,i_to_values,i_to_basis)
end

# Optimization for
#
#  g = lazy_map( linear_combination, cell_to_i_to_val, cell_to_i_to_f)
#  lazy_map(Broadcasting(gradient),g)
#
function lazy_map(
  k::Broadcasting{typeof(∇)}, a::LazyArray{<:Fill{typeof(linear_combination)}})

  i_to_basis = lazy_map(k,a.args[2])
  i_to_values = a.args[1]
  lazy_map(linear_combination,i_to_values,i_to_basis)
end

function lazy_map(
  k::Broadcasting{typeof(∇∇)}, a::LazyArray{<:Fill{typeof(linear_combination)}})

  i_to_basis = lazy_map(k,a.args[2])
  i_to_values = a.args[1]
  lazy_map(linear_combination,i_to_values,i_to_basis)
end

# Optimization for
#
#  g = lazy_map(transpose,cell_to_i_to_f)
#  lazy_map(Broadcasting(gradient),g)
#
function lazy_map(
  k::Broadcasting{typeof(∇)}, a::LazyArray{<:Fill{typeof(transpose)}})

  i_to_basis = lazy_map(k,a.args[1])
  lazy_map( transpose, i_to_basis)
end

function lazy_map(
  k::Broadcasting{typeof(∇∇)}, a::LazyArray{<:Fill{typeof(transpose)}})

  i_to_basis = lazy_map(k,a.args[1])
  lazy_map( transpose, i_to_basis)
end

# Gradient rules
for op in (:+,:-)
  @eval begin

    function lazy_map(
      ::typeof(gradient), a::LazyArray{<:Fill{Operation{typeof($op)}}})

      f = a.args
      g = map(i->lazy_map(gradient,i),f)
      lazy_map(Operation($op),g...)
    end

    function lazy_map(
      ::Broadcasting{typeof(gradient)}, a::LazyArray{<:Fill{Broadcasting{Operation{typeof($op)}}}})

      f = a.args
      g = map(i->lazy_map(Broadcasting(∇),i),f)
      lazy_map(Broadcasting(Operation($op)),g...)
    end

  end
end

for op in (:*,:⋅,:⊙,:⊗)
  @eval begin

    function lazy_map(
      ::typeof(gradient), a::LazyArray{<:Fill{Operation{typeof($op)}}})

      f = a.args
      @notimplementedif length(f) != 2
      g = map(i->lazy_map(gradient,i),f)
      k(F1,F2,G1,G2) = product_rule($op,F1,F2,G1,G2)
      lazy_map(Operation(k),f...,g...)
    end

    function lazy_map(
      ::Broadcasting{typeof(gradient)}, a::LazyArray{<:Fill{Broadcasting{Operation{typeof($op)}}}})

      f = a.args
      @notimplementedif length(f) != 2
      g = map(i->lazy_map(Broadcasting(∇),i),f)
      k(F1,F2,G1,G2) = product_rule($op,F1,F2,G1,G2)
      lazy_map(Broadcasting(Operation(k)),f...,g...)
    end

  end
end

function lazy_map(
  ::Broadcasting{typeof(gradient)}, a::LazyArray{<:Fill{Broadcasting{typeof(∘)}}})

  f = a.args[1]
  g = a.args[2]
  ∇f = lazy_map(Broadcasting(∇),f)
  h = lazy_map(Broadcasting(∘),∇f,g)
  ∇g = lazy_map(Broadcasting(∇),g)

  lazy_map(Broadcasting(Operation(⋅)),∇g,h)
end

function lazy_map(
  ::typeof(gradient), a::LazyArray{<:Fill{typeof(∘)}})

  f = a.args[1]
  g = a.args[2]
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

for op in (:push_∇,:push_∇∇)
  @eval begin
    function lazy_map(
      k::Broadcasting{typeof($op)},
      cell_∇at::LazyArray{<:Fill{typeof(transpose)}},
      cell_map::AbstractArray)

      cell_∇a = cell_∇at.args[1]
      cell_∇b = lazy_map(k,cell_∇a,cell_map)
      cell_∇bt = lazy_map(transpose,cell_∇b)
      cell_∇bt
    end
  end
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

# Memoization

MemoArray(a) = a

#struct MemoArray{T,N,A} <: AbstractArray{T,N}
#  parent::A
#  memo::Dict{Any,Any}
#  function MemoArray(parent::AbstractArray{T,N}) where {T,N}
#    A = typeof(parent)
#    memo = Dict()
#    new{T,N,A}(parent,memo)
#  end
#end
#
## Do not wrap twice.
#MemoArray(parent::MemoArray) = parent
#
#function get_children(n::TreeNode, a::MemoArray)
#   (similar_tree_node(n,a.parent),)
#end
#
#Base.size(a::MemoArray) = size(a.parent)
#Base.axes(a::MemoArray) = axes(a.parent)
#Base.IndexStyle(::Type{MemoArray{T,N,A}}) where {T,N,A} = IndexStyle(A)
#Base.getindex(a::MemoArray,i::Integer) = a.parent[i]
#Base.getindex(a::MemoArray{T,N},i::Vararg{Integer,N}) where {T,N} = a.parent[i...]
#Arrays.array_cache(a::MemoArray) = array_cache(a.parent)
#@inline Arrays.getindex!(cache,a::MemoArray,i::Integer) = getindex!(cache,a.parent,i)
#@inline Arrays.getindex!(cache,a::MemoArray{T,N},i::Vararg{Integer,N}) where {T,N} = getindex!(cache,a.parent,i...)
#Arrays.testitem(a::MemoArray) = testitem(a.parent)
#Arrays.get_array(a::MemoArray) = get_array(a.parent)
#
#function lazy_map(::typeof(evaluate),a::MemoArray,x::AbstractArray{<:Point})
#  key = (:evaluate,objectid(x))
#  if ! haskey(a.memo,key)
#    a.memo[key] = lazy_map(evaluate,a.parent,x)
#  end
#  a.memo[key]
#end
#
#function lazy_map(::typeof(evaluate),a::MemoArray,x::AbstractArray{<:AbstractArray{<:Point}})
#  key = (:evaluate,objectid(x))
#  if ! haskey(a.memo,key)
#    a.memo[key] = lazy_map(evaluate,a.parent,x)
#  end
#  a.memo[key]
#end
#
#function lazy_map(k::typeof(∇),a::MemoArray)
#  lazy_map(Broadcasting(∇),a)
#end
#
#function lazy_map(k::Broadcasting{typeof(∇)},a::MemoArray)
#  key = :gradient
#  if ! haskey(a.memo,key)
#    a.memo[key] = MemoArray(lazy_map(k,a.parent))
#  end
#  a.memo[key]
#end
#
#function lazy_map(
#  k::Broadcasting{typeof(push_∇)},
#  cell_∇a::MemoArray,
#  cell_map::AbstractArray)
#  key = (:push_gradient,objectid(cell_map))
#  if ! haskey(cell_∇a.memo,key)
#    cell_∇a.memo[key] = MemoArray(lazy_map(k,cell_∇a.parent,cell_map))
#  end
#  cell_∇a.memo[key]
#end
#
#function lazy_map(k::typeof(axes),a::MemoArray)
#  lazy_map(k,a.parent)
#end
#
#function lazy_map(k::Reindex{<:MemoArray},::Type{T}, j_to_i::AbstractArray) where T
#  key = (:reindex,objectid(j_to_i))
#  if ! haskey(k.values.memo,key)
#    i_to_v = k.values.parent
#    j_to_v = lazy_map(Reindex(i_to_v),T,j_to_i)
#    k.values.memo[key] = MemoArray(j_to_v)
#  end
#  k.values.memo[key]
#end
#
#function lazy_map(k::PosNegReindex{<:MemoArray,<:MemoArray},::Type{T},i_to_iposneg::AbstractArray) where T
#  values_pos = k.values_pos.parent
#  values_neg = k.values_neg.parent
#  r = lazy_map(PosNegReindex(values_pos,values_neg),i_to_iposneg)
#  MemoArray(r)
#end
#
#function lazy_map(
#  k::typeof(∘), a::MemoArray,b::AbstractArray{<:Field})
#  lazy_map(Broadcasting(∘),a,b)
#end
#
#function lazy_map(
#  k::Broadcasting{typeof(∘)}, a::MemoArray,b::AbstractArray{<:Field})
#  key = (:compose,objectid(b))
#  if ! haskey(a.memo,key)
#    a.memo[key] = MemoArray(lazy_map(k,a.parent,b))
#  end
#  a.memo[key]
#end
#
#function lazy_map(
#  k::typeof(transpose), a::MemoArray)
#  key = :transpose
#  if ! haskey(a.memo,key)
#    a.memo[key] = MemoArray(lazy_map(k,a.parent))
#  end
#  a.memo[key]
#end
