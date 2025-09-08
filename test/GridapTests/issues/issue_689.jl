module Issue689

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.CellData
using Test

function main(transf,T)

  𝒯 = CartesianDiscreteModel((0,1,0,1),(10,10)) |> transf
  labels = get_face_labeling(𝒯)
  add_tag_from_tags!(labels,"surface",6)
  Ω = Interior(𝒯)
  Γ = Boundary(𝒯,tags="surface")

  k = 2
  dΩ = Measure(Ω,2*k)
  dΓ = Measure(Γ,2*k)

  reffe = ReferenceFE(lagrangian,T,k)
  V = FESpace(Ω,reffe)
  Q = FESpace(Γ,reffe)

  u((x,y)) = x^k + y^k
  uh = interpolate(u,V)
  qh = interpolate(u,Q)
  Δ_Γ(f) = x->∇∇(f)(x)[1,1]

  tol = 1.0e-9
  @test sqrt(sum( ∫( abs2( Δ(u) - Δ(uh) )  )dΩ )) < tol
  @test sqrt(sum( ∫( abs2( Δ(u) - Δ(uh) )  )dΓ )) < tol
  @test sqrt(sum( ∫( abs2( Δ_Γ(u) - Δ(qh) )  )dΓ )) < tol

end

main(simplexify,Float64)
main(identity,Float64)

end
