module StateLawsTests

using Gridap.Arrays
using Gridap.Arrays: _split
using Gridap.Fields
import Gridap.Fields: kernel_evaluate
import Gridap.Fields: evaluate_field!
import Gridap.Fields: field_cache
using Gridap.Geometry
using Gridap.Geometry: UnimplementedField
import Gridap.Arrays: kernel_cache
import Gridap.Arrays: apply_kernel!
import Gridap.Arrays: kernel_return_type

import Gridap.Arrays: array_cache
import Gridap.Arrays: getindex!
import Gridap.Fields: evaluate_field_array

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
  k = StateLawKernel(op)
  r = apply_to_field_array(T,k,get_array(a),get_arrays(b...)...)
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

#struct EvalStateLaw{K<:StateLawKernel} <: Kernel
#  k::K
#end
#
#function kernel_evaluate(k::EvalStateLaw,q,a,b...)
#  a_q = evaluate_field_array(a,q)
#  apply(k.k,a_q,b...)
#end

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

using Gridap
using Gridap.Geometry: DiscreteModelMock
using Gridap.FESpaces

#model = DiscreteModelMock()
domain = (0,1,0,1)
n = 3
partition = (n,n)
model = CartesianDiscreteModel(domain,partition)

order = 3
V = TestFESpace(model=model,reffe=:Lagrangian,valuetype=Float64,order=order,conformity=:L2)
U = TrialFESpace(V)

v = get_cell_basis(V)
du = get_cell_basis(V)
uh = FEFunction(U,rand(num_free_dofs(U)))

@statelaw function foo(u,s,r)
  w = u
  w, s, r+1
end

trian = Triangulation(model)
degree = 2*order
quad = CellQuadrature(trian,degree)

q = get_coordinates(quad)
s_q = [i*ones(size(qi)) for (i,qi) in enumerate(q)]
a = ArrayOfEvaluatedFields(s_q,q)
test_array_of_fields(a,q,s_q)

s = QPointCellField(0.0,trian,quad)
r = QPointCellField(0.0,trian,quad)

r_q = evaluate(r,q)
@show r_q === r.array.array

wh = foo(uh,s,r)

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


t_Ω = AffineFETerm((v,u) -> v*u, (v)-> v*r, trian, quad)
op = AffineFEOperator(V,U,t_Ω)

rh = solve(op)

writevtk(trian,"trian",nsubcells=order,cellfields=["r"=>rh])

end # module
