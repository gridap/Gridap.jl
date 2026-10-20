module DiracDeltasTests

using Test
using Gridap.TensorValues
using Gridap.Helpers
using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData

using Gridap.Algebra
using Gridap.Fields
using Gridap.FESpaces
using Gridap.ReferenceFEs

domain = (0,1,0,1)
cells = (30,30)
model = CartesianDiscreteModel(domain,cells)

Ω = Triangulation(model)

vfun(x) = x[1] + 3*x[2]
v = CellField(vfun,Ω)

degree = 2
dirac = DiracDelta{1}(model,degree,tags=5)
@test sum(dirac(v)) ≈ 0.5

#using Gridap.Visualization
#writevtk(get_triangulation(δ),"trian_d")
#writevtk(Ω,"trian")

dirac = DiracDelta{0}(model,tags=2)
@test sum(dirac(v)) ≈ 1

dirac = DiracDelta{0}(model,tags=4)
@test sum(dirac(v)) ≈ 4

@test sum(dirac(3.0)) ≈ 3.0
@test sum(dirac(x->2*x)) ≈ VectorValue(2,2)

dirac = DiracDelta{0}(model,tags=[2,4])
@test sum(dirac(v)) ≈ 5

# Tests for DiracDelta at a generic Point in the domain #

p = Point(0.2,0.3)
dirac = DiracDelta(model,p)
@test sum(dirac(v)) ≈ v(p)

@show dirac = DiracDelta(Ω,p)
@test sum(dirac(v)) ≈ v(p)

pvec = [p,π*p]
dirac = DiracDelta(model,pvec)
@test sum(dirac(v)) ≈ sum(v(pvec))

@show dirac = DiracDelta(Ω,pvec)
@test sum(dirac(v)) ≈ sum(v(pvec))

p = Point(0.2,0.3)
dirac = DiracDelta(model,p)
reffe = ReferenceFE(lagrangian,Float64,3)
V = FESpace(model,reffe,conformity=:L2)
V0 = TestFESpace(model,reffe,conformity=:H1,dirichlet_tags="boundary")
U0 = TrialFESpace(V0)
uh = FEFunction(U0,rand(num_free_dofs(U0)))
vh = FEFunction(V,rand(num_free_dofs(V)))
fe_basis = get_fe_basis(U0)
dg_basis = get_fe_basis(V)

@test sum(dirac(uh)) ≈ uh(p)
@test sum(dirac(vh)) ≈ vh(p)
@test norm(sum(dirac(fe_basis)) .- fe_basis(p)) ≈ 0
@test norm(sum(dirac(dg_basis)) .- dg_basis(p)) ≈ 0

pvec = [p,π*p]
dirac = DiracDelta(model,pvec)

@test sum(dirac(uh)) ≈ sum(map(uh,pvec))
@test sum(dirac(vh)) ≈ sum(map(vh,pvec))
@test norm(sum(dirac(fe_basis)) .- sum(map(fe_basis,pvec))) ≈ 0
@test norm(sum(dirac(dg_basis)) .- sum(map(dg_basis,pvec))) ≈ 0

# comparing exact solution with that of GenericDiracDelta

# 1D Test
model = CartesianDiscreteModel((-1.,1.),(2,))
Ω = Triangulation(model)
dΩ = Measure(Ω,2)

# exact solution
function u(x)
  if x[1] < 0.0
    return 0.0
  else
    return -x[1]
  end
end

ucf = CellField(u,Ω)

order = 1
reffe = ReferenceFE(lagrangian,Float64,order)
V = FESpace(model,reffe,conformity=:L2)
V0 = TestFESpace(model,reffe,conformity=:H1,dirichlet_tags="boundary")
U0 = TrialFESpace(V0,u)

p = Point(0.0)

δ_p = DiracDelta(model,p)
δ_tag = DiracDelta{0}(model,tags=3)

a(u,v) = ∫(∇(u)⋅∇(v))*dΩ
l(v) = ∫(0.0*v)*dΩ + δ_p(1.0*v)

op = AffineFEOperator(a,l,U0,V0)
uh_p = solve(op)

l(v) = ∫(0.0*v)*dΩ + δ_tag(1.0*v)
op = AffineFEOperator(a,l,U0,V0)
uh_tag = solve(op)

e_p = u - uh_p
e = uh_p - uh_tag
err_p = ∫(e_p*e_p)*dΩ
err = ∫(e*e)*dΩ
@test sum(err_p) < eps()
@test sum(err) < eps()

# 2D Test - comparing with tag version and point version

model = CartesianDiscreteModel((-1.,1.,-1.,1.),(3,3))
Ω = Triangulation(model)
dΩ = Measure(Ω,2)

pvec = [Point(-1.0/3,-1.0/3), Point(-1.0/3,1.0/3), Point(1.0/3,1.0/3), Point(1.0/3,-1.0/3)]
δ_p = DiracDelta(model,pvec)

a(u,v) = ∫(∇(u)⋅∇(v))*dΩ
l(v) = δ_p(v)
V0 = TestFESpace(model,reffe,conformity=:H1,dirichlet_tags="boundary")
U0 = TrialFESpace(V0,u)
op = AffineFEOperator(a,l,U0,V0)
uh_p = solve(op)

δ_tag = DiracDelta{0}(model,tags=9)
l(v) = δ_tag(v)
op = AffineFEOperator(a,l,U0,V0)
uh_tag = solve(op)

err = uh_p - uh_tag
@test sum(∫(err*err)*dΩ) < eps()

# testing the faster and direct functional evalulation method

f(x) = sin(norm(x))
fcf = CellField(f,Ω)
@test sum(δ_p(f)) ≈ sum(δ_p(fcf))


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
