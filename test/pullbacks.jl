
using FillArrays
using Gridap
using Gridap.ReferenceFEs, Gridap.FESpaces

using Gridap.ReferenceFEs: Pullback

model = CartesianDiscreteModel((0,2,0,4),(2,2))
Ω = Triangulation(model)
dΩ = Measure(Ω,2)
pts = CellData.get_data(get_cell_points(dΩ))

reffe = RaviartThomasRefFE(Float64,QUAD,1)

V = FESpace(model,reffe)

u(x) = VectorValue(x[1], -x[2])
uh = interpolate(u,V)

φ_phys = get_fe_basis(V).cell_basis
φ_ref  = φ_phys.args[1]
Jt, sign = φ_phys.args[2:end]

σ_phys = get_fe_dof_basis(V).cell_dof
σ_ref  = Fill(get_dof_basis(reffe),num_cells(model))

App = lazy_map(evaluate,σ_phys,φ_phys)[1]
Arr = lazy_map(evaluate,σ_ref,φ_ref)[1]

pf  = ContraVariantPiolaMap()

ipf = inverse_map(pf)
f_ref = lazy_map(ipf,φ_phys,Jt,sign)
Brr = lazy_map(evaluate,σ_ref,f_ref)[1]
Brr == Arr
φ_ref_x = lazy_map(evaluate,φ_ref,pts)[1]
f_ref_x = lazy_map(evaluate,f_ref,pts)[1]
f_ref_x == φ_ref_x

pb  = Pullback(pf)
θ_ref = lazy_map(pb,σ_phys,Jt,sign)

ipb = inverse_map(pb)
