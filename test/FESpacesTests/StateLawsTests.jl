module StateLawsTests

using Gridap.Arrays
using Gridap.Arrays: _split
using Gridap.Fields
import Gridap.Fields: kernel_evaluate
using Gridap.Geometry
using Gridap.Geometry: UnimplementedField
import Gridap.Arrays: kernel_cache
import Gridap.Arrays: apply_kernel!
import Gridap.Arrays: kernel_return_type

"""
"""
function get_arrays(a,b...)
  (get_array(a),get_arrays(b...)...)
end

function get_arrays(a)
  (get_array(a),)
end


"""
"""
@inline function getitems(a::Tuple{Vararg{<:AbstractArray}},i...)
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

function apply_statelaw(op::Function,a::CellField,b...)
  T = UnimplementedField
  k = EvalStateLaw(StateLawKernel(op))
  r = apply(T,k,get_array(a),get_arrays(b...)...)
  similar_object(a,r)
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

struct EvalStateLaw{K<:StateLawKernel} <: Kernel
  k::K
end

function kernel_evaluate(k::EvalStateLaw,q,a,b...)
  a_q = evaluate_field_array(a,q)
  apply(k.k,a_q,b...)
end

using Gridap
using Gridap.Geometry: DiscreteModelMock
using Gridap.FESpaces

model = DiscreteModelMock()

order = 3
V = TestFESpace(model=model,reffe=:Lagrangian,valuetype=Float64,order=order)
U = TrialFESpace(V)

v = get_cell_basis(V)
du = get_cell_basis(V)
uh = FEFunction(U,rand(num_free_dofs(U)))

function foo(u,s,r)
  w = u
  w, s, r+1
end

trian = Triangulation(model)
degree = order
quad = CellQuadrature(trian,degree)

ids = get_cell_id(trian)

q = get_coordinates(quad)
s_q = [i*ones(size(qi)) for (i,qi) in enumerate(q)]
r_q = [zeros(size(qi)) for qi in q]

wh = apply_statelaw(foo,uh,s_q,r_q)

wh_q = evaluate(wh,q)
uh_q = evaluate(uh,q)

#display(uh_q)
#display(wh_q)

function loop(a)
cache = array_cache(a)
for k in 1:length(a)
  ak = getindex!(cache,a,k)
  @show ak
end
end

loop(wh_q)

#display(s_q)
display(r_q)

@show typeof(wh_q)



end # module
