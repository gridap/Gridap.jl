
mass(u,v,dΩ) = ∫(u⋅v)dΩ
laplacian(u,v,dΩ) = ∫(∇(u)⊙∇(v))dΩ

function bm_matrix_assembly(
  D :: Integer,
  n :: Integer,
  reffe :: Tuple,
  qdegree :: Integer,
  biform :: Function
)
  domain = Tuple(repeat((0,1), D)...)
  partition = Tuple(repeat(n, D)...)
  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain, partition))
  Ω = Triangulation(model)
  dΩ = Measure(Ω,qdegree)
  V = TestFESpace(model, reffe)
  a(u,v) = biform(u,v,dΩ)
  A = assemble_matrix(a, V, V)
end
