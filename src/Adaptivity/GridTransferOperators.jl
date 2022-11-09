
abstract type GridTransferOperator <: GridapType end

"""
  ProjectionTransferOperator

  Matrix-like operator providing a transfer between two `FESpaces`
  defined in child-parent adapted grids. 

  The dofs are transfered by iteratively solving a projection problem 

        (uh_to,vh_to)_to = (uh_from,vh_to)_Ω

  where Ω is always chosen as the finest grid (i.e the RHS is
  always integrated without loss of accuracy). 
"""
struct ProjectionTransferOperator{T,A,B,C} <: GridTransferOperator
  from   ::A
  to     ::B
  caches ::C

  function ProjectionTransferOperator(T,from,to,caches)
    A = typeof(from)
    B = typeof(to)
    C = typeof(caches)
    new{T,A,B,C}(from,to,caches)
  end
end

# L2 projection (default for the transfer ops)
Π_l2(u,v,dΩ) = ∫(v⋅u)*dΩ

function ProjectionTransferOperator(from::FESpace,to::FESpace; solver::LinearSolver=BackslashSolver(), Π::Function=Π_l2, qdegree::Int=3)
  @assert isa(from,TrialFESpace)
  @assert isa(to,TrialFESpace)

  Ω_from = get_triangulation(from)
  Ω_to   = get_triangulation(to)
  @assert isa(Ω_from,AdaptedTriangulation) || isa(Ω_to,AdaptedTriangulation)
  @assert is_change_possible(Ω_from,Ω_to)

  # Choose integration space
  Ω  = best_target(Ω_from,Ω_to)
  dΩ = Measure(Ω,qdegree)
  U  = (Ω === Ω_from) ? from : to
  V  = U.space
  vh_to = get_fe_basis(to)
  vh = change_domain(vh_to,Ω,ReferenceDomain())

  # Prepare system
  lhs_mat, lhs_vec = assemble_lhs(Π,Ω_to,to,to.space,qdegree)
  rhs_vec = similar(lhs_vec)
  assem   = SparseMatrixAssembler(to,to.space)

  # Prepare solver
  ss = symbolic_setup(solver,lhs_mat)
  ns = numerical_setup(ss,lhs_mat)

  caches = ns, lhs_vec, rhs_vec, Π, assem, Ω, dΩ, U, V, vh
  return ProjectionTransferOperator(eltype(lhs_mat),from,to,caches)
end

# Solves the problem Π(uh,vh)_to = Π(uh_from,vh)_Ω for all vh in Vh_to
function LinearAlgebra.mul!(y,A::ProjectionTransferOperator,x)
  ns, lhs_vec, rhs_vec, Π, assem, Ω, dΩ, U, V , vh_Ω = A.caches
  Ω_to = get_triangulation(A.to)

  # Bring uh to the integration domain
  uh_from = FEFunction(A.from,x)
  uh_Ω    = change_domain(uh_from,Ω,ReferenceDomain())

  # Assemble rhs vector
  contr   = Π(uh_Ω,vh_Ω,dΩ)
  if Ω !== Ω_to
    contr = merge_contr_cells(contr,Ω,Ω_to)
  end
  vecdata = collect_cell_vector(A.to.space,contr)
  assemble_vector!(rhs_vec,assem,vecdata)
  rhs_vec .-= lhs_vec

  # Solve projection
  solve!(y,ns,rhs_vec)
  return y
end

function Base.size(A::ProjectionTransferOperator)
  (num_free_dofs(A.to),num_free_dofs(A.from))
end

function Base.size(A::ProjectionTransferOperator,i::Int)
  if i == 1
    return num_free_dofs(A.to)
  elseif i == 2
    return num_free_dofs(A.from)
  else
    return nothing
  end
end

function Base.display(op::ProjectionTransferOperator{T,A,B,C}) where {T,A,B,C}
  s = size(op)
  println("$(s[1])x$(s[2])  ProjectionTransferOperator{$(T)}")
end

function assemble_lhs(Π,Ω,Uh,Vh,qdegree)
  dΩ      = Measure(Ω,qdegree)
  uh_dir  = FEFunction(Uh,zero_free_values(Uh),get_dirichlet_dof_values(Uh))
  a(u,v)  = Π(u,v,dΩ)
  b(v)    = a(uh_dir,v)

  lhs_mat = assemble_matrix(a,Uh,Vh)
  lhs_vec = assemble_vector(b,Vh)
  return lhs_mat, lhs_vec
end


# These two functions can actually be rewritten by overloading
# CombineContributionsMap & move_contributions() 
function merge_contr_cells(a::DomainContribution,rtrian::AdaptedTriangulation,ctrian)
  b = DomainContribution()
  for trian in get_domains(a)
    cell_vec = get_contribution(a,trian)
    res = f2c_cell_contrs(rtrian,cell_vec)
    add_contribution!(b,ctrian,res)
  end
  return b
end

function f2c_cell_contrs(trian::AdaptedTriangulation{Dc,Dp},cell_vec) where {Dc,Dp}
  @check num_cells(trian) == length(cell_vec)

  model = get_adapted_model(trian)
  glue  = get_adaptivity_glue(model)
  nC    = num_cells(get_parent(model))
  ccell_to_fcell = glue.c2f_faces_map

  # Map that sums fine contributions for each coarse cell
  return lazy_map((I,V)->sum(V[I]),ccell_to_fcell,Fill(cell_vec,nC))
end
