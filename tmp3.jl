using Gridap
using Gridap.Geometry

function run()

  domain = (0,1,0,1)
  partition = (3,1)
  order = 1
  model_Ω = CartesianDiscreteModel(domain,partition)
  Ω = Triangulation(model_Ω)
  labels = get_face_labeling(model_Ω)
  bgface_to_mask = get_face_mask(labels,[3,4,6],1)
  Γface_to_bgface = findall(bgface_to_mask)
  model_Γ = BoundaryDiscreteModel(Polytope{1},model_Ω,Γface_to_bgface)
  Γ = Triangulation(model_Γ)
  Λb = SkeletonTriangulation(model_Γ)
  dΛb = Measure(Λb,2*order)

  reffe = ReferenceFE(lagrangian,Float64,order)
  V_Ω = TestFESpace(model_Ω,reffe,conformity=:H1)
  V_Γ = TestFESpace(model_Γ,reffe,conformity=:H1)
  U_Ω = TrialFESpace(V_Ω)
  U_Γ = TrialFESpace(V_Γ)
  X = MultiFieldFESpace([U_Ω,U_Γ])
  Y = MultiFieldFESpace([V_Ω,V_Γ])

  # Weak form
  a((ϕ,η),(w,v)) = ∫( v.⁺*η.⁺ )dΛb
  l((w,v)) = 0.0
  op = AffineFEOperator(a,l,X,Y)
  sol = solve(op)

end
