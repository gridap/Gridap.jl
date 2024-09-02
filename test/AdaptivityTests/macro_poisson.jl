using Gridap
using Gridap.CellData, Gridap.Adaptivity
import Gridap.ReferenceFEs: LagrangianRefFE, Quadrature
import Gridap.Adaptivity: RedRefinementRule, MacroReferenceFE, CompositeQuadrature
import Gridap.Adaptivity: get_polytope, num_subcells
import FillArrays: Fill

l2_norm(a) = sum(∫(a⋅a)*dΩ)
l2_error(a,b) = l2_norm(a-b)

order, ncell = 2, 64
bc_tags = Dict(:dirichlet_tags => "boundary")

u((x, y)) = sin(3.2x * (x - y))cos(x + 4.3y) + sin(4.6 * (x + 2y))cos(2.6(y - 2x))
f(x) = -Δ(u)(x)

model = CartesianDiscreteModel((0, 1, 0, 1), (ncell, ncell))
u_reffe = ReferenceFE(lagrangian, Float64, order)
U = TrialFESpace(FESpace(model, u_reffe; bc_tags...), u)

rrule = RefinementRule(QUAD,2)
#rrule = RedRefinementRule(QUAD)
poly = get_polytope(rrule)
reffe = LagrangianRefFE(Float64, poly, 1)
sub_reffes = Fill(reffe, num_subcells(rrule))
v_reffe = MacroReferenceFE(rrule, sub_reffes)
V = FESpace(model, v_reffe; bc_tags...)

Ω = Triangulation(model)
qdegree = 2*order
macro_quad = Quadrature(poly, CompositeQuadrature(), rrule, 2order)
dΩ = Measure(Ω, macro_quad)
dΩ⁺ = Measure(Ω, 10)

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
V = FESpace(model, v_reffe)
vh = interpolate(sol, V)
tprint(vh)
∇vh = gradient(vh)
tprint(∇vh)

U = FESpace(model, u_reffe)
uh = interpolate(sol, U)
tprint(uh)
∇uh = gradient(uh)
tprint(∇uh)


function j(r)
  ∇r = ∇(r)
  ∫(r*r + ∇r⋅∇r)dΩ
end

∇jv = Gridap.gradient(j, vh)
tprint(∇jv)

∇ju = Gridap.gradient(j, uh)
tprint(∇ju)

@benchmark Gridap.gradient($j, $vh)
@benchmark Gridap.gradient($j, $uh)

cell_v = get_cell_dof_values(vh)
cfv = CellField(V,cell_v)
tprint(cfv)

cell_u = get_cell_dof_values(uh)
cfu = CellField(U,cell_u)
tprint(cfu)

fields_v = get_data(cfv)
fields_u = get_data(cfu)
pts = get_cell_points(dΩ.quad)
#pts = get_cell_points(Ω)
pts_data = get_data(pts)
fields_vx = lazy_map(evaluate,fields_v,pts_data)
fields_ux = lazy_map(evaluate,fields_u,pts_data)
@benchmark lazy_map(evaluate,$fields_v,$pts_data)
@benchmark lazy_map(evaluate,$fields_u,$pts_data)
tprint(fields_vx)
tprint(fields_ux)
@which lazy_map(evaluate,fields_v,pts_data)

jvx = get_array(j(cfv))
jvu = get_array(j(cfu))
@benchmark j($cfv)
@benchmark j($cfu)
tprint(jvx)

∇cfv = ∇(cfv)
∇cfu = ∇(cfu)
@benchmark ∇($cfv)
@benchmark ∇($cfu)
tprint(∇cfv)
tprint(∇cfu)


@benchmark ($cfv⋅$cfv)($pts)
@benchmark ($cfu⋅$cfu)($pts)

tprint(∇cfv)
@benchmark ($∇cfv)($pts)
@benchmark ($∇cfu)($pts)

@benchmark lazy_map(evaluate,$∇cfv.cell_field.args[2].args[2],$pts_data)

@benchmark return_cache($∇cfv.cell_field.args[2].args[2].value,$pts_data[1])


_a = ∇cfv.cell_field.args[2].args[2].value
_b = pts_data[1]
caches = map(return_cache,_a.fine_data,_b.fine_data)

@benchmark map(return_cache,_a.fine_data,_b.fine_data)

@benchmark return_cache(_a.fine_data[1],_b.fine_data[1])

@benchmark return_cache(_a.fine_data[1].args[1],_b.fine_data[1])
@benchmark return_cache(_a.fine_data[1].args[2],_b.fine_data[1])

T = eltype(evaluate!(first(caches),first(a.fine_data),first(b.fine_data)))

@which lazy_map(Broadcasting(∇),fields_v)
@which lazy_map(Broadcasting(∇),fields_u)

@benchmark lazy_map(Broadcasting(∇),fields_v)
@benchmark lazy_map(Broadcasting(∇),fields_u)

function jj(r)
  ∇r = ∇(r)
  r*r + ∇r⋅∇r
end

jj(cfv)
jj(cfu)
@benchmark jj($cfv)
@benchmark jj($cfu)


cell_map = get_cell_map(rrule)

cell_map = get_cell_map(get_grid(rrule.ref_grid))



rrule_bis = Gridap.Adaptivity.RefinementRule(QUAD,2)
grid = Gridap.Geometry.compute_reference_grid(QUAD,(2,2))

cell_map = Gridap.Arrays.collect1d(get_cell_map(rrule_bis.ref_grid.grid))
