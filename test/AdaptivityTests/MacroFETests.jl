
using Gridap
using Gridap.Adaptivity, Gridap.Geometry, Gridap.ReferenceFEs, Gridap.Arrays, Gridap.Helpers
using FillArrays

using Gridap.Fields

using Gridap.Adaptivity: num_subcells

Dc = 2
order = 2

model = CartesianDiscreteModel((0,1,0,1),(2,2))
fmodel = refine(model,2)

rrule = Adaptivity.RedRefinementRule(QUAD)
ncells = num_subcells(rrule)

reffe = LagrangianRefFE(Float64,QUAD,order)
sub_reffes = Fill(reffe,ncells)
macro_reffe = Adaptivity.MacroReferenceFE(rrule,sub_reffes)
macro_quad  = Quadrature(QUAD,Adaptivity.CompositeQuadrature(),rrule,2*order)

V = FESpace(model,reffe)
Vm = FESpace(model,macro_reffe)
Vf  = FESpace(fmodel,reffe)
@assert num_free_dofs(Vm) == num_free_dofs(Vf)

u_exact(x) = cos(2π*x[1])*sin(2π*x[2])

# TODO: um is a LinearCombination of FineToCoarseFields. 
# It should be a FineToCoarseField of LinearCombinations...
u = interpolate(u_exact,V)
um = interpolate(u_exact,Vm)
uf = interpolate(u_exact,Vf)

Ω = Triangulation(model)
Ωf = Triangulation(fmodel)

dΩ = Measure(Ω,2*order)
dΩm = Measure(Ω,macro_quad)
dΩf = Measure(Ωf,2*order)
dΩfc = Measure(Ω,Ωf,2*order)

mass(u,v,dΩ) = ∫(u*v)dΩ
M = assemble_matrix((u,v) -> mass(u,v,dΩ),V,V)

bf = assemble_vector(v -> ∫(v*uf)dΩfc,V)
