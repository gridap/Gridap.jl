module DiracDeltasTests

using Test
using Gridap.TensorValues
using Gridap.Helpers
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData

domain = (0,1,0,1)
cells = (30,30)
model = CartesianDiscreteModel(domain,cells)

Ω = Triangulation(model)

vfun(x) = x[1] + 3*x[2]
v = CellField(vfun,Ω)

degree = 2
δ = DiracDelta{1}(model,degree,tags=5)
@test sum(δ(v)) ≈ 0.5

#using Gridap.Visualization
#writevtk(get_triangulation(δ),"trian_d")
#writevtk(Ω,"trian")

δ = DiracDelta{0}(model,tags=2)
@test sum(δ(v)) ≈ 1

δ = DiracDelta{0}(model,tags=4)
@test sum(δ(v)) ≈ 4

@test sum(δ(3.0)) ≈ 3.0
@test sum(δ(x->2*x)) ≈ VectorValue(2,2)

δ = DiracDelta{0}(model,tags=[2,4])
@test sum(δ(v)) ≈ 5

#using Gridap
#
#order = 2
#reffe = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
#V = FESpace(model,reffe,dirichlet_tags=[1,2,5])
#
#dΩ = Measure(Ω,2*order)
#
#δ = DiracDelta{0}(model,tags=4)
#
#const E = 1
#const ν = 0.33
#const λ = (E*ν)/((1+ν)*(1-2*ν))
#const μ = E/(2*(1+ν))
#σ(ε) = λ*tr(ε)*one(ε) + 2*μ*ε
#
#a(u,v) = ∫( ε(v) ⊙ (σ∘ε(u)) )*dΩ
#l(v) = δ( VectorValue(-1,0)⋅v )
#
#op = AffineFEOperator(a,l,V,V)
#uh = solve(op)
#
#writevtk(Ω,"results",cellfields=["uh"=>uh])
#writevtk(δ.Γ,"d",cellfields=["uh"=>uh])

end # module
