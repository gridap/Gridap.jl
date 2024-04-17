
using BenchmarkTools
using SparseArrays, LinearAlgebra
using Gridap, Gridap.FESpaces, Gridap.Algebra

function nzindex(A::SparseArrays.AbstractSparseMatrixCSC, i0::Integer, i1::Integer)
  if !(1 <= i0 <= size(A, 1) && 1 <= i1 <= size(A, 2)); throw(BoundsError()); end
  ptrs = SparseArrays.getcolptr(A)
  r1 = Int(ptrs[i1])
  r2 = Int(ptrs[i1+1]-1)
  (r1 > r2) && return -1
  r1 = searchsortedfirst(rowvals(A), i0, r1, r2, Base.Order.Forward)
  ((r1 > r2) || (rowvals(A)[r1] != i0)) ? 0 : r1
end

function precompute_nzindex(A,a::Algebra.AllocationCOO)
  precompute_nzindex(A,a.I,a.J)
end

function precompute_nzindex(A,I,J)
  K = zeros(Int32,length(I))
  for (p,(i,j)) in enumerate(zip(I,J))
      if i < 1 || j < 1
          continue
      end
      K[p] = nzindex(A,i,j)
  end
  K
end

function sparse_matrix!(A,V,K;reset=true)
  if reset
      LinearAlgebra.fillstored!(A,0)
  end
  A_nz = nonzeros(A)
  for (k,v) in zip(K,V)
      if k < 1
          continue
      end
      A_nz[k] += v
  end
  A
end

function Algebra.CounterCOO(::Algebra.SparseMatrixBuilder{T},axes) where T
  Algebra.CounterCOO{T}(axes)
end

function assemble_coo(a,matdata;reuse=false)
  m1 = Algebra.CounterCOO(get_matrix_builder(a),(get_rows(a),get_cols(a)))
  symbolic_loop_matrix!(m1,a,matdata)
  m2 = nz_allocation(m1)
  numeric_loop_matrix!(m2,a,matdata)
  m3 = create_from_nz(m2)
  if reuse
    return m3, (m2,precompute_nzindex(m3,m2))
  end
  m3
end

function assemble_coo!(mat,a,matdata,cache)
  m2, K = cache
  m2.counter.nnz = 0
  numeric_loop_matrix!(m2,a,matdata)
  sparse_matrix!(mat,m2.V,K,reset=true)
end

model = CartesianDiscreteModel((0,1,0,1),(100,100))

reffe = ReferenceFE(lagrangian, Float64, 1)
V = TestFESpace(model, reffe)

Ω = Triangulation(model)
dΩ = Measure(Ω, 2)
a(u,v) = ∫(∇(u)⋅∇(v))dΩ

u = get_trial_fe_basis(V)
v = get_fe_basis(V)
matdata = collect_cell_matrix(V,V,a(u,v))

assem = SparseMatrixAssembler(V,V)

A1 = assemble_matrix(assem,matdata)
A2, cache = assemble_coo(assem,matdata;reuse=true)
assemble_coo!(A2,assem,matdata,cache)

@benchmark assemble_matrix($assem,$matdata)
@benchmark assemble_coo($assem,$matdata)
@benchmark assemble_matrix!($A1,$assem,$matdata)
@benchmark assemble_coo!($A2,$assem,$matdata,$cache)
