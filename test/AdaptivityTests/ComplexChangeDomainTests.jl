module ComplexChangeDomainTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using FillArrays

norm2(u,Ω::Triangulation) = sum(∫(u⋅u)*Measure(Ω,qorder))
norm2(u,dΩ::Measure) = sum(∫(u⋅u)*dΩ)

order = 1
qorder = 2*order+1
sol(x) = sum(x)

D = 2
domain = Tuple(repeat([0,1],D))
cmodel = CartesianDiscreteModel(domain,Tuple(fill(2,D)))
fmodel = refine(cmodel,2)

Ωc = Triangulation(cmodel)
Ωf = Triangulation(fmodel)

reffe = ReferenceFE(lagrangian,Float64,order)
Vc = TestFESpace(cmodel,reffe)
Uc = TrialFESpace(Vc)
Vf = TestFESpace(fmodel,reffe)
Uf = TrialFESpace(Vf)

"""
  BodyFittedTriangulation Views
"""

Ωc_view = view(Ωc,[1,2])
Ωf_view = view(Ωf,[1,2,3,4,5,6,7,8])

v_Ωf = interpolate(sol,Vf)
v_Ωf_view = change_domain(v_Ωf,Ωf_view,ReferenceDomain())

v_Ωc_view_view = change_domain(v_Ωf_view,Ωc_view,ReferenceDomain())

v_Ωc = interpolate(sol,Vc)
v_Ωc_view = change_domain(v_Ωc,Ωc_view,ReferenceDomain())

v_Ωf_view_view = change_domain(v_Ωc_view,Ωf_view,ReferenceDomain())

@test norm2(v_Ωc_view_view,Ωc_view) ≈ norm2(v_Ωc_view,Ωc_view)
@test norm2(v_Ωf_view_view,Ωf_view) ≈ norm2(v_Ωf_view,Ωf_view)

# Same but automatic
@test norm2(v_Ωc_view,Ωf_view) ≈ norm2(v_Ωf_view,Ωf_view)
@test norm2(v_Ωf_view,Ωc_view) ≈ norm2(v_Ωc_view,Ωc_view)

"""
  BodyFitted -> Boundary
"""
Γc = Boundary(cmodel)
Γf = Boundary(fmodel)

v_Ωc = interpolate(sol,Vc)
v_Ωf = interpolate(sol,Vf)

v_Γf = change_domain(v_Ωf,Γf,ReferenceDomain())
v_Γc = change_domain(v_Ωc,Γc,ReferenceDomain())

@test norm2(v_Ωf,Γc) ≈ norm2(v_Γc,Γc)
@test norm2(v_Ωc,Γf) ≈ norm2(v_Γf,Γf)

"""
  BodyFitted -> Skeleton

  Not currently working, we need to treat SkeletonPairs separately....
"""

Λc = Skeleton(cmodel)
Λf = Skeleton(fmodel)

for sign in [:plus,:minus]
  v_Ωc = getproperty(interpolate(sol,Vc),sign)
  v_Ωf = getproperty(interpolate(sol,Vf),sign)

  v_Λf = change_domain(v_Ωf,Λf,ReferenceDomain())
  v_Λc = change_domain(v_Ωc,Λc,ReferenceDomain())

  @test norm2(v_Ωf,Λc) ≈ norm2(v_Λc,Λc)
  @test norm2(v_Ωc,Λf) ≈ norm2(v_Λf,Λf)
end

end