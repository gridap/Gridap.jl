
using FillArrays
using Gridap
using Gridap.ReferenceFEs, Gridap.FESpaces, Gridap.CellData
using Gridap.Fields

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

# Inverse Pushforward
ipf = inverse_map(pf)
f_ref = lazy_map(ipf,φ_phys,Jt,sign)
Brr = lazy_map(evaluate,σ_ref,f_ref)[1]
Brr == Arr
φ_ref_x = lazy_map(evaluate,φ_ref,pts)[1]
f_ref_x = lazy_map(evaluate,f_ref,pts)[1]
f_ref_x ≈ φ_ref_x

# Pullback
pb  = Pullback(pf)
θ_ref = lazy_map(pb,σ_phys,Jt,sign)
θ_ref_x = lazy_map(evaluate,θ_ref,φ_ref)
σ_ref_x = lazy_map(evaluate,σ_ref,φ_ref)
θ_ref_x ≈ σ_ref_x

# Inverse Pullback
ipb = inverse_map(pb)
θ_phys = lazy_map(ipb,σ_ref,Jt,sign)
θ_phys_x = lazy_map(evaluate,θ_phys,φ_phys)
σ_phys_x = lazy_map(evaluate,σ_phys,φ_phys)
θ_phys_x ≈ σ_phys_x
