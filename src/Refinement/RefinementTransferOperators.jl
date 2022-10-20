
"""
  RefinementTransferOperator

  Matrix-like operator providing a transfer between two `FESpaces`
  defined in child-parent refined grids. 

  The dofs are transfered by iteratively solving a projection problem 

        (uh_to,vh_to)_to = (uh_from,vh_to)_Ω

  where Ω is always chosen as the finest grid (i.e the RHS is
  always integrated without loss of accuracy). 
"""
struct RefinementTransferOperator{T,A,B,C} <: AbstractMatrix{T}
  from   ::A
  to     ::B
  caches ::C

  function RefinementTransferOperator(T,from,to,caches)
    A = typeof(from)
    B = typeof(to)
    C = typeof(caches)
    new{T,A,B,C}(from,to,caches)
  end
end

# L2 projection (default for the transfer ops)
Π_l2(u,v,dΩ) = ∫(v⋅u)*dΩ

function RefinementTransferOperator(from::FESpace,to::FESpace; Π=Π_l2, qdegree=3)
  @assert isa(from,TrialFESpace)
  @assert isa(to,TrialFESpace)

  Ω_from = get_triangulation(from)
  Ω_to   = get_triangulation(to)
  @assert isa(Ω_from,RefinedTriangulation) || isa(Ω_to,RefinedTriangulation)
  @assert is_change_possible(Ω_from,Ω_to)

  # Choose integration space
  Ω  = best_target(Ω_from,Ω_to)
  dΩ = Measure(Ω,qdegree)
  U  = (Ω === Ω_from) ? from : to
  V  = U.space
  vh_to = get_fe_basis(to)
  vh = change_domain(vh_to,Ω)

  # Prepare system
  sysmat, sysvec = assemble_lhs(Π,Ω_to,to,to.space,qdegree)
  assem  = SparseMatrixAssembler(to,to.space)

  cache = sysmat, sysvec, Π, assem, Ω, dΩ, U, V, vh
  return RefinementTransferOperator(eltype(sysmat),from,to,cache)
end

# Solves the problem Π(uh,vh)_to = Π(uh_from,vh)_Ω for all vh in Vh_to
function LinearAlgebra.mul!(y,A::RefinementTransferOperator,x)
  sysmat, sysvec, Π, assem, Ω, dΩ, U, V , vh_Ω = A.caches
  Ω_to = get_triangulation(A.to)

  # Bring uh to the integration domain
  uh_from = FEFunction(A.from,x)
  uh_Ω    = change_domain(uh_from,Ω)

  # Assemble rhs vector
  contr   = Π(uh_Ω,vh_Ω,dΩ)
  if Ω !== Ω_to
    contr = merge_contr_cells(contr,Ω,Ω_to)
  end
  vecdata = collect_cell_vector(A.to.space,contr)
  assemble_vector_add!(sysvec,assem,vecdata)

  # Solve projection
  IterativeSolvers.cg!(y,sysmat,sysvec)
  return y
end

function Base.size(A::RefinementTransferOperator)
  (num_free_dofs(A.to),num_free_dofs(A.from))
end

function Base.size(A::RefinementTransferOperator,i::Int)
  if i == 1
    return num_free_dofs(A.to)
  elseif i == 2
    return num_free_dofs(A.from)
  else
    return nothing
  end
end

function Base.display(op::RefinementTransferOperator{T,A,B,C}) where {T,A,B,C}
  s = size(op)
  println("$(s[1])x$(s[2])  RefinementTransferOperator{$(T)}")
end

function assemble_lhs(Π,Ω,Uh,Vh,qdegree)
  dΩ = Measure(Ω,qdegree)
  uh_dir = FEFunction(Uh,zero_free_values(Uh),get_dirichlet_dof_values(Uh))
  a(u,v) = Π(u,v,dΩ)
  b(v)   = a(uh_dir,v)

  sysmat = assemble_matrix(a,Uh,Vh)
  sysvec = assemble_vector(b,Vh)
  return sysmat, -sysvec
end

function merge_contr_cells(a::DomainContribution,rtrian::RefinedTriangulation,ctrian)
  b = DomainContribution()
  for trian in get_domains(a)
    cell_vec = get_contribution(a,trian)
    res = f2c_cell_contrs(rtrian,cell_vec)
    add_contribution!(b,ctrian,res)
  end
  return b
end

function f2c_cell_contrs(trian::RefinedTriangulation{Dc,Dp},cell_vec) where {Dc,Dp}
  @check num_cells(trian) == length(cell_vec)

  model = get_refined_model(trian)
  glue  = get_glue(model)
  nF    = num_cells(trian)
  nC    = num_cells(get_parent(model))
  nChildren = 4 # Num fine cells per coarse cell

  # Invert fcell_to_ccell
  fcell_to_ccell = glue.f2c_faces_map[Dc+1]
  ccell_to_fcell = [fill(-1,nChildren) for i in 1:nC]
  cidx = fill(1,nC)
  for iF in 1:nF
    iC = fcell_to_ccell[iF]
    ccell_to_fcell[iC][cidx[iC]] = iF
    cidx[iC] += 1
  end

  # Map that sums fine contributions for each coarse cell
  return lazy_map((I,V)->sum(V[I]),ccell_to_fcell,Fill(cell_vec,nC))
end


"""
"""
struct RefinementTransferMap{A<:FESpace,B<:FESpace,C<:RefinementTransferOperator} <: Map
  from ::A
  to   ::B
  op   ::C
end

function RefinementTransferMap(from::FESpace,to::FESpace; Π=Π_l2, qdegree=3)
  op = RefinementTransferOperator(from,to;Π=Π,qdegree=qdegree)
  return RefinementTransferMap(from,to,op)
end

function Arrays.return_cache(m::RefinementTransferMap,uh::FEFunction)
  y = zeros(size(m.op,1))
  return y
end

function Arrays.evaluate!(cache,m::RefinementTransferMap,uh::FEFunction)
  @check get_fe_space(uh) === m.from
  y = cache
  x = get_free_dof_values(uh)
  mul!(y,m.op,x)
  return FEFunction(m.to,y)
end
