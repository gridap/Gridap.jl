
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

###########
mutable struct CachedInserterCOO{Tv,Ti,Ti2}
  nnz :: Int
  mat :: SparseMatrixCSC{Tv,Ti}
  K :: Vector{Ti2}

  function CachedInserterCOO(mat::SparseMatrixCSC{Tv,Ti}, K::Vector{Ti2}) where {Tv,Ti,Ti2}
    new{Tv,Ti,Ti2}(0,mat,K)
  end
end

Algebra.LoopStyle(::Type{<:CachedInserterCOO}) = Loop()

@noinline function Algebra.add_entry!(::typeof(+),a::CachedInserterCOO,v,i,j)
  a.nnz = a.nnz + 1
  k = a.K[a.nnz]
  a.mat.nzval[k] += v
  nothing
end

function assemble_coo_opt(a,matdata;reuse=false)
  m1 = Algebra.CounterCOO(get_matrix_builder(a),(get_rows(a),get_cols(a)))
  symbolic_loop_matrix!(m1,a,matdata)
  m2 = nz_allocation(m1)
  numeric_loop_matrix!(m2,a,matdata)
  m3 = create_from_nz(m2)
  if reuse
    K = precompute_nzindex(m3,m2)
    return m3, K
  end
  m3
end

function assemble_coo_opt!(mat,a,matdata,K)
  LinearAlgebra.fillstored!(mat,zero(eltype(mat)))
  m = CachedInserterCOO(mat,K)
  numeric_loop_matrix!(m,a,matdata)
  mat
end

model = CartesianDiscreteModel((0,1,0,1),(100,100))

order = 1
reffe = ReferenceFE(lagrangian, Float64, order)
V = TestFESpace(model, reffe)

qdegree = 2*order
Ω  = Triangulation(model)
dΩ = Measure(Ω, qdegree)
Γ  = Boundary(model)
dΓ = Measure(Γ, qdegree)
a(u,v) = ∫(∇(u)⋅∇(v))dΩ + ∫(u*v)dΓ - ∫((∇⋅u) ⋅ (∇⋅v))dΩ

u = get_trial_fe_basis(V)
v = get_fe_basis(V)
matdata = collect_cell_matrix(V,V,a(u,v))

assem = SparseMatrixAssembler(V,V)

A1 = assemble_matrix(assem,matdata)
A2, cache2 = assemble_coo(assem,matdata;reuse=true);
assemble_coo!(A2,assem,matdata,cache2)
A3, cache3 = assemble_coo_opt(assem,matdata;reuse=true);
assemble_coo_opt!(A3,assem,matdata,cache3)
A1 ≈ A2 ≈ A3

@benchmark assemble_matrix($assem,$matdata)
@benchmark assemble_coo($assem,$matdata)
@benchmark assemble_coo_opt($assem,$matdata)
@benchmark assemble_matrix!($A1,$assem,$matdata)
@benchmark assemble_coo!($A2,$assem,$matdata,$cache2)
@benchmark assemble_coo_opt!($A3,$assem,$matdata,$cache3)