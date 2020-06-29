const CONT = 0
const DISC = 1

function CDLagrangianRefFE(::Type{T},p::Polytope{D},order::Int,cont) where {T,D}
  orders = tfill(order,Val{D}())
  CDLagrangianRefFE(T,p,orders,cont)
end

function CDLagrangianRefFE(::Type{T},p::Polytope{D},orders,cont) where {T,D}

  lrfe = LagrangianRefFE(T,p,orders)
  grfe = lrfe.data

  face_own_nodes = _compute_cd_face_own_nodes(p,orders,cont)

  face_own_dofs = _generate_face_own_dofs(face_own_nodes, grfe.dofs.node_and_comp_to_dof)

  reffaces = compute_lagrangian_reffaces(T,p,orders)

  _reffaces = vcat(reffaces...)

  ndofs = length(grfe.dofs.dof_to_node)

  face_own_dofs_permutations = _trivial_face_own_dofs_permutations(face_own_dofs)

  GenericRefFE(
      grfe.ndofs,
      grfe.polytope,
      grfe.prebasis,
      grfe.dofs,
      face_own_dofs,
      face_own_dofs_permutations,
      grfe.face_dofs, 
      grfe.shapefuns)

end

function _compute_cd_face_own_nodes(p::ExtrusionPolytope{D},orders::NTuple{D,<:Int},cont::NTuple{D,<:Int}) where D
  nodes_ijk = CartesianIndices(Tuple(orders.+1))
  nodes_i = LinearIndices(nodes_ijk)
  face_owned_nodes = Vector{Int64}[]
  for nf in p.dface.nfaces
    anc = nf.anchor
    ext = nf.extrusion
    fns = UnitRange{Int64}[]
    for (i,(e,a,c,o)) in enumerate(zip(ext,anc,cont,orders))
      if e==0 && c == DISC
        push!(fns,1:-1)
        break
      elseif e==0 && c == CONT
        push!(fns,a*o+1:a*o+1)
      elseif e==1 && c == DISC
        push!(fns,1:o+1)
      elseif e==1 && c == CONT
        push!(fns,2:o)
      end
    end
    n_ijk = CartesianIndices(Tuple(fns))
    if length(n_ijk) == 0
      n_i = Int64[]
    else
      n_i = Int64[nodes_i[n_ijk]...]
    end
    push!(face_owned_nodes,n_i)
  end
  face_owned_nodes
end
