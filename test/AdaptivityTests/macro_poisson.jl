using Gridap
import Gridap.ReferenceFEs: LagrangianRefFE, Quadrature
import Gridap.Adaptivity: RedRefinementRule, MacroReferenceFE, CompositeQuadrature
import Gridap.Adaptivity: get_polytope, num_subcells
import FillArrays: Fill


order, ncell = 2, 64
bc_tags = Dict(:dirichlet_tags => "boundary")

u((x, y)) = sin(3.2x * (x - y))cos(x + 4.3y) + sin(4.6 * (x + 2y))cos(2.6(y - 2x))
f(x) = -Δ(u)(x)

model = CartesianDiscreteModel((0, 1, 0, 1), (ncell, ncell))
u_reffe = ReferenceFE(lagrangian, Float64, order)
U = TrialFESpace(FESpace(model, u_reffe; bc_tags...), u)

rrule = RedRefinementRule(QUAD)
poly = get_polytope(rrule)
reffe = LagrangianRefFE(Float64, poly, 1)
sub_reffes = Fill(reffe, num_subcells(rrule))
v_reffe = MacroReferenceFE(rrule, sub_reffes)
V = FESpace(model, v_reffe; bc_tags...)

Ω = Triangulation(model)
macro_quad = Quadrature(poly, CompositeQuadrature(), rrule, 2order)
dΩ = Measure(Ω, macro_quad)
dΩ⁺ = Measure(Ω, 10)

a(u, v) = ∫(∇(v) ⋅ ∇(u))dΩ
l(v) = ∫(f * v)dΩ
op = AffineFEOperator(a, l, U, V)
eh = solve(op) - u
fem_l2err = sqrt(∑(∫(eh * eh)dΩ))    # 0.5745, too large!

a1(u, v) = ∫(u * v)dΩ
l1(v) = ∫(u * v)dΩ
op1 = AffineFEOperator(a1, l1, U, V)
eh1 = solve(op1) - u
fem_l2err1 = sqrt(∑(∫(eh1 * eh1)dΩ))  # 1.52e-5, reasonable

############################################################################################

l2_norm(a) = sum(∫(a⋅a)*dΩ)
l2_error(a,b) = l2_norm(a-b)

sol(x) = sum(x)
V = FESpace(model, v_reffe)

uh = interpolate(sol, V)
l2_error(uh, sol)

∇uh = gradient(uh)
l2_error(∇uh, ∇(sol))

Ωrr = Triangulation(rrule.ref_grid)

cmaps = get_cell_map(rrule)
cjacs = lazy_map(∇,cmaps)
cinvjacs = lazy_map(inv,cjacs)

pts = Gridap.CellData.get_data(get_cell_points(Ωrr))

lazy_map(evaluate,cinvjacs,pts)
