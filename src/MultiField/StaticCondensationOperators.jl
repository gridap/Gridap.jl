
struct StaticCondensationOperator <: FEOperator
  full_space :: FESpace
  eliminated_space :: FESpace
  retained_space :: FESpace
  patch_assem :: Assembler
  eliminated_assem :: Assembler
  retained_assem :: Assembler
  full_matvecs
  sc_op
end

function StaticCondensationOperator(
  full_space :: FESpace,
  eliminated_space :: FESpace,
  retained_space :: FESpace,
  patch_assem,
  full_matvecs,
)
  n_elim = num_fields(eliminated_space)
  n_ret = num_fields(retained_space)
  @assert isa(MultiFieldStyle(full_space),BlockMultiFieldStyle{2,(n_elim,n_ret)})

  eliminated_assem = SparseMatrixAssembler(eliminated_space,eliminated_space)
  retained_assem = SparseMatrixAssembler(retained_space,retained_space)

  Asc, bsc = statically_condensed_assembly(retained_assem,patch_assem,full_matvecs)
  sc_op = AffineFEOperator(retained_space,retained_space,Asc,bsc)

  return StaticCondensationOperator(
    full_space,eliminated_space,retained_space,
    patch_assem,eliminated_assem,retained_assem,
    full_matvecs,sc_op
  )
end

function StaticCondensationOperator(
  full_space::FESpace, args...; kwargs...
)
  @check isa(MultiFieldStyle(full_space),BlockMultiFieldStyle{2})
  eliminated_space, retained_space = blocks(full_space)
  return StaticCondensationOperator(
    full_space, eliminated_space, retained_space, args...; kwargs...
  )
end

function StaticCondensationOperator(
  ptopo,
  full_space :: FESpace,
  eliminated_space :: FESpace,
  retained_space :: FESpace,
  biform :: Function,
  liform :: Function
)
  patch_assem = FESpaces.PatchAssembler(ptopo,full_space,full_space)
  full_matvecs = assemble_matrix_and_vector(biform,liform,patch_assem,full_space,full_space,zero(full_space))
  return StaticCondensationOperator(
    full_space,eliminated_space,retained_space,
    patch_assem,full_matvecs
  )
end

function StaticCondensationOperator(
  ptopo, full_space::FESpace, args...; kwargs...
)
  @check isa(MultiFieldStyle(full_space),BlockMultiFieldStyle{2})
  eliminated_space, retained_space = blocks(full_space)
  return StaticCondensationOperator(
    ptopo, full_space, eliminated_space, retained_space, args...; kwargs...
  )
end

function Algebra.solve(op::StaticCondensationOperator)
  solver = LinearFESolver()
  solve(solver,op)
end

function Algebra.solve!(uh,solver::LinearFESolver,op::StaticCondensationOperator,::Nothing)
  x_eliminated, x_retained = blocks(get_free_dof_values(uh))
  uh_eliminated = FEFunction(op.eliminated_space,x_eliminated)
  uh_retained = FEFunction(op.retained_space,x_retained)
  cache = solve!(uh_retained,solver,op.sc_op,nothing)
  backward_static_condensation!(uh_eliminated,op,uh_retained)
  return uh, cache
end

function Algebra.solve!(uh,solver::LinearFESolver,op::StaticCondensationOperator,cache)
  x_eliminated, x_retained = blocks(get_free_dof_values(uh))
  uh_eliminated = FEFunction(op.eliminated_space,x_eliminated)
  uh_retained = FEFunction(op.retained_space,x_retained)
  cache = solve!(uh_retained,solver,op.sc_op,cache)
  backward_static_condensation!(uh_eliminated,op,uh_retained)
  return uh, cache
end

function statically_condensed_assembly(retained_assem,patch_assem,full_matvecs)
  sc_matvecs = lazy_map(StaticCondensationMap(),full_matvecs)
  rows = patch_assem.strategy.array.array[2,2].patch_rows
  cols = patch_assem.strategy.array.array[2,2].patch_cols
  data = (([sc_matvecs,],[rows,],[cols,]), ([],[],[]), ([],[]))
  assemble_matrix_and_vector(retained_assem,data)
end

function backward_static_condensation!(x_eliminated, eliminated_assem,patch_assem,full_matvecs,x_retained)
  rows_elim = patch_assem.strategy.array.array[1,1].patch_rows
  rows_ret = patch_assem.strategy.array.array[2,2].patch_rows

  patch_x_ret = lazy_map(Broadcasting(Reindex(x_retained)),rows_ret)
  patch_x_elim = lazy_map(BackwardStaticCondensationMap(),full_matvecs,patch_x_ret)

  vecdata = ([patch_x_elim,],[rows_elim,])
  assemble_vector!(x_eliminated,eliminated_assem,vecdata)
end

function backward_static_condensation!(
  x_eliminated::AbstractVector,op::StaticCondensationOperator,x_retained::AbstractVector
)
  backward_static_condensation!(x_eliminated,op.eliminated_assem,op.patch_assem,op.full_matvecs,x_retained)
end

function backward_static_condensation!(
  xh_eliminated,op::StaticCondensationOperator,xh_retained
)
  x_eliminated = get_free_dof_values(xh_eliminated)
  x_retained = get_free_dof_values(xh_retained)
  backward_static_condensation!(x_eliminated,op,x_retained)
  return xh_eliminated
end

function backward_static_condensation(op::StaticCondensationOperator,x_retained::AbstractVector)
  x_eliminated = zero_free_values(op.eliminated_space)
  backward_static_condensation!(x_eliminated,op,x_retained)
  return x_eliminated
end

function backward_static_condensation(op::StaticCondensationOperator,xh_retained)
  xh_eliminated = zero(op.eliminated_space)
  backward_static_condensation!(xh_eliminated,op,xh_retained)
  return xh_eliminated
end

FESpaces.get_trial(op::StaticCondensationOperator) = op.full_space
FESpaces.get_test(op::StaticCondensationOperator) = op.full_space

struct StaticCondensationMap{A} <: Map
  pivot :: A
end

StaticCondensationMap() = StaticCondensationMap(RowMaximum())

function Arrays.evaluate!(cache,k::StaticCondensationMap, matvec::Tuple)
  mat, vec = matvec
  evaluate!(cache,k,mat,vec)
end

function Arrays.evaluate!(cache,k::StaticCondensationMap, mat, vec)
  @check size(mat) == (2,2)
  @check size(vec) == (2,)

  Kii, Kbi, Kib, Kbb = get_array(mat)
  bi, bb = get_array(vec)

  f = lu!(Kii,k.pivot;check=false)
  @check issuccess(f) "Factorization failed"
  ldiv!(f,bi)
  ldiv!(f,Kib)

  mul!(bb,Kbi,bi,-1,1)
  mul!(Kbb,Kbi,Kib,-1,1)

  return Kbb, bb
end

struct BackwardStaticCondensationMap{A} <: Map
  pivot :: A
end

BackwardStaticCondensationMap() = BackwardStaticCondensationMap(RowMaximum())

function Arrays.evaluate!(cache,k::BackwardStaticCondensationMap, matvec::Tuple, xb)
  mat, vec = matvec
  evaluate!(cache,k,mat,vec,xb)
end

function Arrays.evaluate!(cache,k::BackwardStaticCondensationMap, mat, vec, xb)
  @check size(mat) == (2,2)
  @check size(vec) == (2,)

  Kii, Kbi, Kib, Kbb = get_array(mat)
  bi, bb = get_array(vec)

  f = lu!(Kii, k.pivot; check=false)
  @check issuccess(f) "Factorization failed"

  # Reconstruct interior solution
  mul!(bi, Kib, xb, -1, 1)  # bi = bi - Kib * xb
  ldiv!(f, bi)              # bi = Kii^{-1} * (bi - Kib * xb)

  return bi
end
