const CONT = 0
const DISC = 1

struct CDConformity{D} <: Conformity
  cont::NTuple{D,Int}
end

function Conformity(reffe::GenericLagrangianRefFE{<:CDConformity},sym::Symbol)
  if sym == :L2
    L2Conformity()
  else
    @unreachable """\n
    It is not possible to use conformity = $sym on a LagrangianRefFE with CD conformity.

    Either use conformity = :L2 or pass an instance of CDConformity.
    """
  end
end

function get_face_own_nodes(reffe::GenericLagrangianRefFE{<:CDConformity},conf::CDConformity)
  _cd_get_face_own_nodes(reffe,conf)
end

function get_face_own_nodes(reffe::GenericLagrangianRefFE{GradConformity},conf::CDConformity)
  _cd_get_face_own_nodes(reffe,conf)
end

# Constructors

function _CDLagrangianRefFE(::Type{T},p::ExtrusionPolytope{D},order::Int,cont) where {T,D}
  orders = tfill(order,Val{D}())
  _CDLagrangianRefFE(T,p,orders,cont)
end

function _CDLagrangianRefFE(::Type{T},p::ExtrusionPolytope{D},orders,cont) where {T,D}
  cond(c,o) = ( o > 0 || c == DISC )
  @check all((cond(cont[k],orders[k]) for k in 1:length(orders)))
  _cd_lagrangian_ref_fe(T,p,orders,cont)
end

function _cd_lagrangian_ref_fe(::Type{T},p::ExtrusionPolytope{D},orders,cont) where {T,D}

  @check isa(p,ExtrusionPolytope)

  prebasis = compute_monomial_basis(T,p,orders)

  nodes, face_own_nodes = cd_compute_nodes(p,orders)
  dofs = LagrangianDofBasis(T,nodes)

  nnodes = length(dofs.nodes)
  ndofs = length(dofs.dof_to_node)

  face_own_nodes = _compute_cd_face_own_nodes(p,orders,cont)
  face_nodes = _compute_face_nodes(p,face_own_nodes)

  face_own_dofs = _generate_face_own_dofs(face_own_nodes, dofs.node_and_comp_to_dof)
  face_dofs = _compute_face_nodes(p,face_own_dofs)

  data = nothing

  conf = CDConformity(Tuple(cont))

  reffe = GenericRefFE{typeof(conf)}(
      ndofs,
      p,
      prebasis,
      dofs,
      conf,
      data,
      face_dofs)

  GenericLagrangianRefFE(reffe,face_nodes)
end

function _cd_get_face_own_nodes(reffe,conf::CDConformity)
  p = get_polytope(reffe)
  orders = get_orders(get_prebasis(reffe))
  cont = conf.cont
  cond(c,o) = ( o > 0 || c == DISC )
  @assert all((cond(cont[k],orders[k]) for k in 1:length(orders)))
  p = get_polytope(reffe)
  dofs = get_dof_basis(reffe)
  @assert is_n_cube(p)
  face_own_nodes = _compute_cd_face_own_nodes(p,orders,cont)
  face_own_nodes
end

function _cd_get_face_own_dofs(reffe,conf::CDConformity)
  p = get_polytope(reffe)
  orders = get_orders(get_prebasis(reffe))
  cont = conf.cont
  cond(c,o) = ( o > 0 || c == DISC )
  @assert all((cond(cont[k],orders[k]) for k in 1:length(orders)))
  p = get_polytope(reffe)
  dofs = get_dof_basis(reffe)
  @assert is_n_cube(p)
  face_own_nodes = _compute_cd_face_own_nodes(p,orders,cont)
  face_own_dofs = _generate_face_own_dofs(face_own_nodes, dofs.node_and_comp_to_dof)
  face_own_dofs
end

function cd_compute_nodes(p::Polytope{D},orders) where D
  active_faces = _active_faces(p,orders)
  nodes = Point{D,Float64}[]
  facenodes = [Int[] for i in 1:num_faces(p)]
  _cd_compute_high_order_nodes_dim_0!(nodes,facenodes,p,active_faces)
  for d in 1:(num_dims(p)-1)
    _cd_compute_high_order_nodes_dim_d!(nodes,facenodes,p,orders,Val{d}(),active_faces)
  end
  if active_faces[end]
    _cd_compute_high_order_nodes_dim_D!(nodes,facenodes,p,orders)
  end

  nodes, facenodes = _move_nodes(p,orders,active_faces,nodes,facenodes)

  (nodes, facenodes)
end

function _cd_compute_high_order_nodes_dim_0!(nodes,facenodes,p,active_faces)
  x = get_vertex_coordinates(p)
  k = 1
  for vertex in 1:num_vertices(p)
    if active_faces[vertex]
      push!(nodes,x[vertex])
      push!(facenodes[vertex],k)
      k += 1
    end
  end
end

@noinline function _cd_compute_high_order_nodes_dim_d!(nodes,facenodes,p,orders,::Val{d},active_faces) where d
  x = get_vertex_coordinates(p)
  offset = get_offset(p,d)
  k = length(nodes)+1
  for iface in 1:num_faces(p,d)
    if active_faces[offset+iface]
      face = Polytope{d}(p,iface)
      face_ref_x = get_vertex_coordinates(face)
      face_prebasis = MonomialBasis(Float64,face,1)
      change = inv(evaluate(face_prebasis,face_ref_x))
      face_shapefuns = linear_combination(change,face_prebasis)
      face_vertex_ids = get_faces(p,d,0)[iface]
      face_x = x[face_vertex_ids]
      face_orders = compute_face_orders(p,face,iface,orders)
      face_interior_nodes = compute_own_nodes(face,face_orders)
      face_high_x = evaluate(face_shapefuns,face_interior_nodes)*face_x
      for xi in 1:length(face_high_x)
        push!(nodes,face_high_x[xi])
        push!(facenodes[iface+offset],k)
        k += 1
      end
    end
  end
end

function _cd_compute_high_order_nodes_dim_D!(nodes,facenodes,p,orders)
  k = length(nodes)+1
  p_high_x = compute_own_nodes(p,orders)
  for xi in 1:length(p_high_x)
    push!(nodes,p_high_x[xi])
    push!(facenodes[end],k)
    k += 1
  end
end

function _compute_cd_face_own_nodes(p::ExtrusionPolytope{D},orders::NTuple{D,<:Int},cont::NTuple{D,<:Int}) where D
  nodes, _ = cd_compute_nodes(p,orders)
  _ord = convert(NTuple{D,Float64},orders)
  _nodes = map(i -> i.data.*_ord,nodes)
  _nodes = map(i -> convert(NTuple{D,Int},i),_nodes)
  perm = Dict(_nodes[i] => i for i = 1:length(_nodes))
  face_owned_nodes = Vector{Int}[]
  for nf in p.dface.nfaces
    anc = nf.anchor
    ext = nf.extrusion
    fns = UnitRange{Int}[]
    for i in 1:length(orders)
      e,a,c,o = ext[i],anc[i],cont[i],orders[i]
      if e==0 && c == DISC
        push!(fns,0:-1)
        break
      elseif e==0 && c == CONT
        push!(fns,a*o:a*o)
      elseif e==1 && c == DISC
        push!(fns,0:o)
      elseif e==1 && c == CONT
        push!(fns,1:o-1)
      end
    end
    n_ijk = map(i->i.I,collect(CartesianIndices(Tuple(fns))))
    if length(n_ijk) == 0
      n_i = Int[]
    else
      n_i = sort(reshape([perm[i] for i in n_ijk],length(n_ijk)))
    end
    push!(face_owned_nodes,n_i)
  end
  face_owned_nodes
end

function _compute_face_nodes(p::ExtrusionPolytope,f_o)
  nfs = p.dface.nf_nfs
  f_n = [Int[] for i in 1:num_faces(p)]
  for i in 1:length(nfs)
    for nf in nfs[i]
      push!(f_n[i],f_o[nf]...)
    end
  end
  f_n
end

function _active_faces(p::Polytope,orders)
  nfs = p.dface.nfaces
  is_valid = Bool[]
  cond = (o,e,a) -> (o > 0 || (e == 0 && a == 0))
  for nf in nfs
    isv = all((cond(orders[k],nf.extrusion[k],nf.anchor[k]) for k in 1:length(orders)))
    push!(is_valid,isv)
  end
  is_valid
end

function _move_nodes(p::ExtrusionPolytope,orders,act_f,nodes,f_own)
  if any(map(i->i==0,orders))
    tr = map(i -> i == 0 ? 1 : 0,orders)
    nfs = p.dface.nfaces
    con(ext,anc) = f -> ( f.extrusion == ext && f.anchor == anc )
    perm = zeros(Int,length(nfs))
    for i in findall(act_f)
      nf  =nfs[i]
      ext = VectorValue(map(+,nf.extrusion,tr))
      anc = nf.anchor
      r = findall(con(ext,anc),nfs)
      @assert (length(r)<=1)
      perm[i] = r[1]
    end

    facenodes = [Int[] for i in 1:num_faces(p)]
    for f in findall(act_f)
      facenodes[perm[f]] = f_own[f]
    end

    tr = VectorValue(map(i -> i == 0 ? 0.5 : 0.0,orders))
    nodes = nodes.+tr

    return nodes, facenodes

  else
    return nodes, f_own
  end
end
