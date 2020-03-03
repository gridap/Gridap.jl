
"""
"""
function apply_statelaw(op::Function,a::CellField,b...)
  T = UnimplementedField
  k = StateLawKernel(op)
  r = apply_to_field_array(T,k,get_array(a),get_arrays(b...)...)
  similar_object(a,r)
end

"""
"""
macro statelaw(fundef)
  s = "The @statelaw macro is only allowed in function definitions"
  @assert isa(fundef,Expr) s
  @assert fundef.head in (:(=), :function) s
  funname = fundef.args[1].args[1]
  nargs = length(fundef.args[1].args)-1
  @assert nargs >= 2 "The @statelaw macro is only allowed in functions with at leats two arguments"

  y = fundef.args[1].args[2]
  x =  fundef.args[1].args[3:end]
  a = fundef.args[1].args[2:end]
  q = quote
    function $(funname)($(y)::CellField,$(x...))
      apply_statelaw($(funname),$(a...))
    end
    $(fundef)
  end
  :($(esc(q)))
end

struct StateLawKernel{F<:Function} <: Kernel
  op::F
end

function kernel_cache(k::StateLawKernel,a::AbstractVector,b::AbstractVector...)
  v = zeros(eltype(a),size(a))
  CachedArray(v)
end

function kernel_return_type(k::StateLawKernel,a::AbstractVector,b::AbstractVector...)
  v = zeros(eltype(a),size(a))
  typeof(v)
end

function apply_kernel!(cache,k::StateLawKernel,a::AbstractVector,b::AbstractVector...)
  Q = length(a)
  setsize!(cache,size(a))
  v = cache.array
  for q in 1:Q
    aq = a[q]
    bq = getitems(b,q)
    r = k.op(aq,bq...)
    vq, states = _split(r...)
    _update_states!(b,q,states,Val{length(states)}())
    v[q] = vq
  end
  v
end

function apply_kernel_for_cache!(cache,k::StateLawKernel,a::AbstractVector,b::AbstractVector...)
  Q = length(a)
  setsize!(cache,size(a))
  v = cache.array
  for q in 1:Q
    aq = a[q]
    bq = getitems(b,q)
    r = k.op(aq,bq...)
    vq, states = _split(r...)
    v[q] = vq
  end
  v
end

@inline function _update_states!(b,q,states,::Val{i}) where i
  _update_state!(b,q,states,Val{i}())
  _update_states!(b,q,states,Val{i-1}())
  nothing
end

@inline function _update_states!(b,q,states,::Val{0})
  nothing
end

@inline function _update_state!(b,q,states,::Val{i}) where i
  m = length(b)
  n = length(states)
  o = m-n
  b[i+o][q] = states[i]
  nothing
end

