
using Test
using Gridap
using Gridap.Adaptivity, Gridap.Geometry, Gridap.ReferenceFEs, Gridap.Arrays, Gridap.Helpers
using Gridap.CellData, Gridap.Fields, Gridap.FESpaces
using FillArrays

using Gridap.Fields

using Gridap.Adaptivity: num_subcells

Dc = 2
order = 2

model = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(2,2)))
fmodel = refine(model)

rrule = Adaptivity.RedRefinementRule(QUAD)
ncells = num_subcells(rrule)

reffe = LagrangianRefFE(Float64,QUAD,order)
sub_reffes = Fill(reffe,ncells)
macro_reffe = Adaptivity.MacroReferenceFE(rrule,sub_reffes)
macro_quad  = Quadrature(QUAD,Adaptivity.CompositeQuadrature(),rrule,2*order)

Ω = Triangulation(model)
Ωf = Triangulation(fmodel)

dΩ = Measure(Ω,2*order)
dΩm = Measure(Ω,macro_quad)
dΩf = Measure(Ωf,2*order)
dΩfc = Measure(Ω,Ωf,2*order)

pt = get_cell_points(dΩ)
ptm = get_cell_points(dΩm)
ptf = get_cell_points(dΩf)

V = FESpace(model,reffe)
Vm = FESpace(model,macro_reffe)
Vf  = FESpace(fmodel,reffe)
@assert num_free_dofs(Vm) == num_free_dofs(Vf)

u = get_trial_fe_basis(V)
um = get_trial_fe_basis(Vm)
v = get_fe_basis(V)
vm = get_fe_basis(Vm)

cf1 = vm⋅u
cf1(pt)
cf1(ptm)

cf2 = um⋅v
cf2(pt)
cf2(ptm)

∇v = ∇(vm)(ptm)

# um and uf should be the same
u_exact(x) = cos(2π*x[1])*sin(2π*x[2])
um = interpolate(u_exact,Vm)
uf = interpolate(u_exact,Vf)

M = assemble_matrix((u,v) -> ∫(u*v)dΩ,V,V)
bf = assemble_vector(v -> ∫(v*uf)dΩfc,V)
bm = assemble_vector(v -> ∫(v*um)dΩm,V)
@test norm(bf-bm) < 1.e-10
xf = M\bf
xm = M\bm
@test norm(xf-xm) < 1.e-10

