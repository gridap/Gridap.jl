module bm_assembly

using PkgBenchmark, BenchmarkTools
using Gridap
using Gridap.Geometry

mass(u,v,dΩ) = ∫(u⋅v)dΩ
laplacian(u,v,dΩ) = ∫(∇(u)⊙∇(v))dΩ
graddiv(u,v,dΩ) = ∫((∇⋅u)⋅(∇⋅v))dΩ

function driver(
  Ω       :: Triangulation,
  reffe   :: Tuple,
  qdegree :: Integer,
  biform  :: Function
)
  model = get_background_model(Ω)
  dΩ = Measure(Ω,qdegree)
  V = TestFESpace(model, reffe)
  a(u,v) = biform(u,v,dΩ)
  A = assemble_matrix(a, V, V)
end

const SUITE = BenchmarkGroup()

for (D,n) in [(2,10),(3,6)]
  domain = Tuple([repeat([0,1], D)...])
  partition = Tuple(fill(n, D))
  model = UnstructuredDiscreteModel(CartesianDiscreteModel(domain, partition))
  trian_cases = [
    ("bulk",Triangulation(model)),
    ("view",Triangulation(model,collect(1:div(n^D,2)))), # About half of the cells
  ]
  for (trian_name, trian) in trian_cases
    for order in [1,2,3]
      basis_cases = [
        ("lagrangian",lagrangian,Float64,order),
        ("vector_lagragian",lagrangian,VectorValue{D,Float64},order),
        ("raviart_thomas",raviart_thomas,Float64,order-1),
      ]
      for (basis_name,basis,T,degree) in basis_cases
        biform_cases = [
          ("mass",mass,2*order),
          ("laplacian",laplacian,2*(order-1)),
          ("graddiv",graddiv,2*(order-1)),
        ]
        for (biform_name,biform,qdegree) in biform_cases
          reffe = ReferenceFE(basis, T, degree)
          name = "assembly_$(D)D_$(trian_name)_$(basis_name)_$(biform_name)_$(order)"
          SUITE[name] = @benchmarkable driver(
            $(trian),$(reffe),$(qdegree),$(biform)
          )
        end
      end
    end
  end
end

end # module