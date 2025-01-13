
using Gridap

using Gridap.Geometry
using Gridap.FESpaces

model = CartesianDiscreteModel((0,1,0,1),(4,4))
Ω = Triangulation(model)

order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V = FESpace(model,reffe)

dΩ = Measure(Ω,2*order)
mass(u,v) = ∫(u*v)dΩ



