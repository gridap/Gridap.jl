using BenchmarkTools
using Gridap

include("drivers.jl")

const SUITE = BenchmarkGroup()

ncells = 40
for D in [2,3]
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
      ]
      for (biform_name,biform,qdegree) in biform_cases
        reffe = ReferenceFE(basis, T, degree)
        name = "assembly_$(D)D_$(basis_name)_$(biform_name)_$(order)"
        SUITE[name] = @benchmarkable bm_matrix_assembly(D,ncells,reffe,qdegree,biform)
      end
    end
  end
end
