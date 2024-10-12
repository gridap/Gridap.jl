
using SparseArrays, LinearAlgebra
using Gridap, Gridap.Algebra, Gridap.FESpaces, Gridap.Helpers

struct MinReuse end

struct ReusableCounterCSC{Tv,Ti}
  tv::Type{Tv}
  nrows::Int
  ncols::Int
  colnnzmax::Vector{Ti}
end

Algebra.LoopStyle(::Type{<:ReusableCounterCSC}) = Loop()

@inline function Algebra.add_entry!(::Function,a::ReusableCounterCSC{Tv,Ti},v,i,j) where {Tv,Ti}
  a.colnnzmax[j] += Ti(1)
  nothing
end

mutable struct ReusableInserterCSC{Tv,Ti}
  nrows::Int
  ncols::Int
  colptr::Vector{Ti}
  colnnz::Vector{Ti}
  rowval::Vector{Ti}
  nzval::Vector{Tv}
  nnz::Ti
  nnzptr::Vector{Ti}
end

Algebra.LoopStyle(::Type{<:ReusableInserterCSC}) = Loop()

@noinline function Algebra.add_entry!(::typeof(+),a::ReusableInserterCSC{Tv,Ti},v,i,j) where {Tv,Ti}
  pini = Int(a.colptr[j])
  pend = pini + Int(a.colnnz[j]) - 1
  p = searchsortedfirst(a.rowval,i,pini,pend,Base.Order.Forward)
  if (p>pend)
    # add new entry
    a.colnnz[j] += 1
    a.rowval[p] = i
    a.nzval[p] = v
  elseif a.rowval[p] != i
    # shift one forward from p to pend
    @check  pend+1 < Int(a.colptr[j+1])
    for k in pend:-1:p
      o = k + 1
      a.rowval[o] = a.rowval[k]
      a.nzval[o] = a.nzval[k]
    end
    # add new entry
    a.colnnz[j] += 1
    a.rowval[p] = i
    a.nzval[p] = v
  else
    # update existing entry
    a.nzval[p] += v 
  end
  # Save the position of the entry for reuse
  a.nnz += 1
  a.nnzptr[a.nnz] = p
  nothing
end

############################################################################################

function Algebra.nz_counter(
  ::SparseMatrixBuilder{SparseMatrixCSC{Tv,Ti},<:MinReuse},
  axes
) where {Tv,Ti}
  nrows = length(axes[1])
  ncols = length(axes[2])
  colnnzmax = zeros(Ti,ncols)
  ReusableCounterCSC(Tv,nrows,ncols,colnnzmax)
end

function Algebra.nz_allocation(a::ReusableCounterCSC{Tv,Ti}) where {Tv,Ti}
  colptr = Vector{Ti}(undef,a.ncols+1)
  @inbounds for i in 1:a.ncols
    colptr[i+1] = a.colnnzmax[i]
  end
  length_to_ptrs!(colptr)
  ndata = colptr[end] - one(Ti)
  rowval = Vector{Ti}(undef,ndata)
  nzval = zeros(Tv,ndata)
  nnzmax = sum(a.colnnzmax)
  nnzptr = fill(zero(Ti),nnzmax)
  colnnz = a.colnnzmax
  fill!(colnnz,zero(Ti))
  ReusableInserterCSC(a.nrows,a.ncols,colptr,colnnz,rowval,nzval,0,nnzptr)
end

mutable struct InserterCache{Ti}
  nnz :: Ti
  nnzptr :: Vector{Ti}
end

function Algebra.create_from_nz(a::ReusableInserterCSC{Tv,Ti}) where {Tv,Ti}
  k = 1
  nnzperm = zeros(Ti,a.nnz)
  for j in 1:a.ncols
    pini = Int(a.colptr[j])
    pend = pini + Int(a.colnnz[j]) - 1
    for p in pini:pend
      a.nzval[k] = a.nzval[p]
      a.rowval[k] = a.rowval[p]
      nnzperm[p] = k
      k += 1
    end
  end
  @inbounds for j in 1:a.ncols
    a.colptr[j+1] = a.colnnz[j]
  end
  length_to_ptrs!(a.colptr)
  nnz = a.colptr[end]-1
  resize!(a.rowval,nnz)
  resize!(a.nzval,nnz)
  nnzptr = collect(lazy_map(Reindex(nnzperm),a.nnzptr))
  SparseMatrixCSC(a.nrows,a.ncols,a.colptr,a.rowval,a.nzval), InserterCache(a.nnz,nnzptr)
end

struct CachedInserterCSC{Tv,Ti}
  mat::SparseMatrixCSC{Tv,Ti}
  cache::InserterCache{Ti}
end

Algebra.LoopStyle(::Type{<:CachedInserterCSC}) = Loop()

@noinline function Algebra.add_entry!(::typeof(+),a::CachedInserterCSC,v,i,j)
  a.cache.nnz = a.cache.nnz + 1
  k = a.cache.nnzptr[a.cache.nnz]
  a.mat.nzval[k] += v
  nothing
end

function FESpaces.assemble_matrix!(mat,a::SparseMatrixAssembler,matdata,cache)
  LinearAlgebra.fillstored!(mat,zero(eltype(mat)))
  cache.nnz = 0
  m = CachedInserterCSC(mat,cache)
  numeric_loop_matrix!(m,a,matdata)
  mat
end

############################################################################################

model = CartesianDiscreteModel((0,1,0,1),(10,10))

reffe = ReferenceFE(lagrangian, Float64, 1)
V = TestFESpace(model, reffe)

Ω = Triangulation(model)
dΩ = Measure(Ω, 2)
a(u,v) = ∫(∇(u)⋅∇(v))dΩ

u = get_trial_fe_basis(V)
v = get_fe_basis(V)
matdata = collect_cell_matrix(V,V,a(u,v))

assem = SparseMatrixAssembler(SparseMatrixBuilder(SparseMatrixCSC{Float64,Int},MinReuse()),V,V)
A1, cache = assemble_matrix(assem,matdata)

A2 = deepcopy(A1)
assemble_matrix!(A2,assem,matdata,cache)

A1 ≈ A2
