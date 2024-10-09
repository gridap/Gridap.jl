using Gridap
using Gridap.CellData, Gridap.Adaptivity
import Gridap.ReferenceFEs: LagrangianRefFE, Quadrature
import Gridap.Adaptivity: RedRefinementRule, MacroReferenceFE, CompositeQuadrature
import Gridap.Adaptivity: get_polytope, num_subcells
import FillArrays: Fill
using Gridap.FESpaces

l2_norm(a) = sum(∫(a⋅a)*dΩ)
l2_error(a,b) = l2_norm(a-b)

order, ncell = 2, 100
bc_tags = Dict(:dirichlet_tags => "boundary")

u((x, y)) = sin(3.2x * (x - y))cos(x + 4.3y) + sin(4.6 * (x + 2y))cos(2.6(y - 2x))
f(x) = -Δ(u)(x)

model = CartesianDiscreteModel((0, 1, 0, 1), (ncell, ncell))
u_reffe = ReferenceFE(lagrangian, Float64, order)
U = TrialFESpace(FESpace(model, u_reffe; bc_tags...), u)

rrule = RefinementRule(QUAD,order)
poly = get_polytope(rrule)
reffe = LagrangianRefFE(Float64, poly, 1)
sub_reffes = Fill(reffe, num_subcells(rrule))
v_reffe = MacroReferenceFE(rrule, sub_reffes)
V = FESpace(model, v_reffe; bc_tags...)

Ω = Triangulation(model)
qdegree = 2*order+2
macro_quad = Quadrature(poly, CompositeQuadrature(), rrule, 2)
dΩ = Measure(Ω, macro_quad)
dΩ⁺ = Measure(Ω, qdegree)
quad = Quadrature(QUAD, qdegree)

a(u, v) = ∫(∇(v) ⋅ ∇(u))dΩ
l(v) = ∫(f * v)dΩ
op = AffineFEOperator(a, l, U, V)
eh = solve(op) - u
fem_l2err = sqrt(∑(∫(eh * eh)dΩ))

a1(u, v) = ∫(u * v)dΩ
l1(v) = ∫(u * v)dΩ
op1 = AffineFEOperator(a1, l1, U, V)
eh1 = solve(op1) - u
fem_l2err1 = sqrt(∑(∫(eh1 * eh1)dΩ))

############################################################################################
import Gridap.CellData: get_data
using BenchmarkTools
using Gridap.Fields

tprint(x::AbstractArray) = print_op_tree(x)
tprint(x::CellField) = print_op_tree(get_data(x))
tprint(x::DomainContribution) = print_op_tree(get_array(x))

sol(x) = sum(x)
V = FESpace(model, v_reffe; dirichlet_tags="boundary")
vh = interpolate(sol, V)
tprint(vh)
∇vh = gradient(vh)
tprint(∇vh)

U = FESpace(model, u_reffe; dirichlet_tags="boundary")
uh = interpolate(sol, U)
tprint(uh)
∇uh = gradient(uh)
tprint(∇uh)

function j(r)
  ∇r = ∇(r)
  ∫(r*r + ∇r⋅∇r)dΩ
end

function j_ref(r)
  ∇r = ∇(r)
  ∫(r*r + ∇r⋅∇r)dΩ⁺
end

∇jv = Gridap.gradient(j, vh)
tprint(∇jv)

∇ju = Gridap.gradient(j_ref, uh)
tprint(∇ju)

@benchmark Gridap.gradient($j, $vh)
@benchmark Gridap.gradient($j_ref, $vh)
@benchmark Gridap.gradient($j, $uh)
@benchmark Gridap.gradient($j_ref, $uh)

v_vec = assemble_vector(Gridap.gradient(j, vh), V)
u_vec = assemble_vector(Gridap.gradient(j_ref, uh), U)
@benchmark assemble_vector(Gridap.gradient($j, $vh), $V)
@benchmark assemble_vector(Gridap.gradient($j_ref, $uh), $U)

using Gridap.Arrays
using ForwardDiff
cell_v = get_cell_dof_values(vh)
dummy_forwarddiff_tag = ()->()
i_to_cfg = lazy_map(Arrays.ConfigMap(ForwardDiff.gradient,dummy_forwarddiff_tag),cell_v)
i_to_xdual = lazy_map(Arrays.DualizeMap(ForwardDiff.gradient,dummy_forwarddiff_tag),cell_v)

cfv = CellField(V, i_to_xdual)
@benchmark CellField($V, $i_to_xdual)
@benchmark j($vh)

function jj(r)
  ∇r = ∇(r)
  return r*r + ∇r⋅∇r
end

jj_v = jj(vh)
@benchmark jj($vh)

cfa = vh*vh
@benchmark $vh⋅$vh

∇vh = ∇(vh)
@benchmark ∇($vh)

cfb = ∇vh⋅∇vh
@benchmark $∇vh ⋅ $∇vh
@benchmark $cfa + $cfb

using Gridap.Arrays
args = (cfa,cfb)
x = CellData._get_cell_points(args...)
ax = map(i->i(x),args)
axi = map(first,ax)
r = Fields.BroadcastingFieldOpMap(+)(axi...)

x = zeros(Point{2,Float64},4)
args_data = map(CellData.get_data,args)
args_vals = map(a -> return_value(first(a),x),args_data)
r = Fields.BroadcastingFieldOpMap(+)(args_vals...)

f1 = first(first(args_data))
@which return_value(f1,x)


function g1(cfa,cfb)
  args = (cfa,cfb)
  x = CellData._get_cell_points(args...)
  ax = map(i->i(x),args)
  axi = map(first,ax)
  r = Fields.BroadcastingFieldOpMap(+)(axi...)
end
function g2(cfa,cfb)
  args = (cfa,cfb)
  x = zeros(Point{2,Float64},1)
  args_data = map(CellData.get_data,args)
  args_vals = map(a -> return_value(first(a),x),args_data)
  r = Fields.BroadcastingFieldOpMap(+)(args_vals...)
end
function g3(cfa,cfb)
  args = (cfa,cfb)
  x = CellData.get_data(CellData._get_cell_points(args...))
  data = map(CellData.get_data,args)
  fields = lazy_map(Broadcasting(+),data...)
  vals = lazy_map(evaluate,fields,x)
  r = testitem(vals)
end

@benchmark g1($cfa,$cfb)
@benchmark g2($cfa,$cfb)
@benchmark g3($cfa,$cfb)

tprint(cfv)

cell_u = get_cell_dof_values(uh)
cfu = CellField(U, cell_u)
tprint(cfu)

∇jv_arr = get_array(∇jv)#.args[1]
∇ju_arr = get_array(∇ju)#.args[1]
tprint(∇jv_arr)
tprint(∇ju_arr)

∇jv_arr_cache = array_cache(∇jv_arr)
getindex!(∇jv_arr_cache, ∇jv_arr, 1)
∇ju_arr_cache = array_cache(∇ju_arr)
getindex!(∇ju_arr_cache, ∇ju_arr, 1)

@benchmark getindex!($∇jv_arr_cache, $∇jv_arr, 100)
@benchmark getindex!($∇ju_arr_cache, $∇ju_arr, 100)

function lazy_collect(cache,arr)
  n = length(arr)
  for i in 1:n
    getindex!(cache, arr, i)
  end
end

@benchmark lazy_collect($∇jv_arr_cache, $∇jv_arr)
@benchmark lazy_collect($∇ju_arr_cache, $∇ju_arr)

@benchmark assemble_vector($∇jv, $V)
@benchmark assemble_vector($∇ju, $U)

vecdata_v = collect_cell_vector(V, ∇jv)
vecdata_u = collect_cell_vector(U, ∇ju)
@benchmark collect_cell_vector($V, $∇jv)
@benchmark collect_cell_vector($U, $∇ju)

assem_u = SparseMatrixAssembler(U,U)
assem_v = SparseMatrixAssembler(V,V)
vec_u = zeros(num_free_dofs(U))
vec_v = zeros(num_free_dofs(V))

numeric_loop_vector!(vec_v,assem_v,vecdata_v)
numeric_loop_vector!(vec_u,assem_u,vecdata_u)
@benchmark numeric_loop_vector!($vec_v, $assem_v, $vecdata_v)
@benchmark numeric_loop_vector!($vec_u, $assem_u, $vecdata_u)

using Gridap.Algebra
function numeric_loop_vector_caches(b,a,vecdata)
  strategy = FESpaces.get_assembly_strategy(a)
  (cellvec, _cellids) = first(zip(vecdata...))
  cellids = FESpaces.map_cell_rows(strategy,_cellids)
  rows_cache = array_cache(cellids)
  vals_cache = array_cache(cellvec)
  vals1 = getindex!(vals_cache,cellvec,1)
  rows1 = getindex!(rows_cache,cellids,1)
  add! = FESpaces.AddEntriesMap(+)
  add_cache = return_cache(add!,b,vals1,rows1)
  return add_cache, vals_cache, rows_cache
end

caches_v = numeric_loop_vector_caches(vec_v,assem_v,vecdata_v)
caches_u = numeric_loop_vector_caches(vec_u,assem_u,vecdata_u)

@benchmark numeric_loop_vector_caches($vec_v, $assem_v, $vecdata_v)
@benchmark numeric_loop_vector_caches($vec_u, $assem_u, $vecdata_u)

cellvec_v, cellids_v = first(zip(vecdata_v...))
cellvec_u, cellids_u = first(zip(vecdata_u...))

FESpaces._numeric_loop_vector!(vec_v,caches_v,cellvec_v,cellids_v)
FESpaces._numeric_loop_vector!(vec_u,caches_u,cellvec_u,cellids_u)
@benchmark FESpaces._numeric_loop_vector!($vec_v, $caches_v, $cellvec_v, $cellids_v)
@benchmark FESpaces._numeric_loop_vector!($vec_u, $caches_u, $cellvec_u, $cellids_u)

function _numeric_loop_vector_bis!(vec,caches,cell_vals,cell_rows)
  add_cache, vals_cache, rows_cache = caches
  @assert length(cell_vals) == length(cell_rows)
  add! = FESpaces.AddEntriesMap(+)
  for cell in 1:length(cell_rows)
    rows = getindex!(rows_cache,cell_rows,cell)
    #vals = getindex!(vals_cache,cell_vals,cell)
    #evaluate!(add_cache,add!,vec,vals,rows)
  end
end

_numeric_loop_vector_bis!(vec_v,caches_v,cellvec_v,cellids_v)
_numeric_loop_vector_bis!(vec_u,caches_u,cellvec_u,cellids_u)
@benchmark _numeric_loop_vector_bis!($vec_v, $caches_v, $cellvec_v, $cellids_v)
@benchmark _numeric_loop_vector_bis!($vec_u, $caches_u, $cellvec_u, $cellids_u)

tprint(cellvec_v)
tprint(cellvec_u)

