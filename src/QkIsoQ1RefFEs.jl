struct FineToCoarseBasis{A<:AbstractArray{<:AbstractArray{<:Field}},
                         B<:Gridap.Adaptivity.RefinementRule} <: AbstractVector{Field}
  fine_fields    :: A
  rrule          :: B
  function FineToCoarseBasis(fine_fields,rrule::Gridap.Adaptivity.RefinementRule)
    Gridap.Helpers.@check length(fine_fields)   == Gridap.Adaptivity.num_subcells(rrule)
    A = typeof(fine_fields)
    B = typeof(rrule)
    new{A,B}(fine_fields,rrule)
  end
end

function Base.length(f::FineToCoarseBasis)
  Gridap.Geometry.num_nodes(f.rrule.ref_grid)
end
function Base.size(f::FineToCoarseBasis)
  (Gridap.Geometry.num_nodes(f.rrule.ref_grid),)
end
# Dirty hack
function Base.getindex(f::FineToCoarseBasis,i)
  fields = f.fine_fields
  points = get_node_coordinates(f.rrule.ref_grid)
  fi = fields[1]
  yi_type  = Fields.return_type(fi,points[1])
  ConstantField(zero(eltype(yi_type)))
end

# Methods required by RefFE re. the prebasis
function Gridap.ReferenceFEs.get_orders(a::FineToCoarseBasis)
  D=num_cell_dims(a.rrule.ref_grid)
  nodes=num_vertices(a.rrule.ref_grid)
  orders=zeros(Int,D)
  order=Int64(nodes^(1.0/D))-1
  orders.=order
  Tuple(orders)
end

function Gridap.ReferenceFEs.get_order(a::FineToCoarseBasis)
  maximum(Gridap.ReferenceFEs.get_orders(a))
end


function Gridap.Geometry.return_cache(a::FineToCoarseBasis,x::AbstractArray{<:Point})
  fields, x2cell = a.fine_fields, xi->Gridap.Adaptivity.x_to_cell(a.rrule,xi)
  cmaps = Gridap.Adaptivity.get_inverse_cell_map(a.rrule)

  xi_cache = array_cache(x)
  fi_cache = array_cache(fields)
  mi_cache = array_cache(cmaps)

  xi = getindex!(xi_cache,x,1)
  child_id = x2cell(xi)
  mi = getindex!(mi_cache,cmaps,child_id)
  fi = getindex!(fi_cache,fields,child_id)

  zi_cache = Fields.return_cache(mi,xi)
  zi = evaluate!(zi_cache,mi,xi)

  yi_type  = Fields.return_type(fi,zi)
  yi_cache = Fields.return_cache(fi,zi)
  y_cache  = Gridap.Arrays.CachedArray(eltype(yi_type),2)

  return fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache
end

function Gridap.Geometry.evaluate!(cache,
                            a::FineToCoarseBasis{<:AbstractArray{<:AbstractArray{<:Field}}},
                            x::AbstractArray{<:Point})
  fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache = cache
  fields, x_to_cell = a.fine_fields, xi->Gridap.Adaptivity.x_to_cell(a.rrule,xi)
  cmaps = Gridap.Adaptivity.get_inverse_cell_map(a.rrule)

  Gridap.Arrays.setsize!(y_cache, (length(x),Gridap.Geometry.num_nodes(a.rrule.ref_grid)))

  y_cache.array .= zero(eltype(y_cache.array))

  cell_node_ids = Gridap.Geometry.get_cell_node_ids(a.rrule.ref_grid)
  for i in eachindex(x)
     xi = getindex!(xi_cache,x,i)
     child_id = x_to_cell(xi)
     fi = getindex!(fi_cache,fields,child_id)
     mi = getindex!(mi_cache,cmaps,child_id)
     zi = Fields.evaluate!(zi_cache,mi,xi)
     y_cache.array[i,cell_node_ids[child_id]] = Fields.evaluate!(yi_cache,fi,zi)
  end
  return y_cache.array
end

function QkIsoQ1(::Type{T},D::Integer,order) where {T}
    orders=Tuple(fill(order,D))
    
    MH=CartesianDiscreteModel(Tuple(repeat([0,1],D)), Tuple(fill(1,D)))
    Mh=Gridap.Adaptivity.refine(MH,orders)
    ref_fe=ReferenceFE(lagrangian,T,1)
    Vl=FESpace(Mh,ref_fe)
    vl=get_fe_basis(Vl)

    rrule=Mh.glue.refinement_rules[1]
    
    prebasis = FineToCoarseBasis(Gridap.CellData.get_data(vl),rrule)

    if D==1
        p=SEGMENT
    elseif D==2
        p=QUAD
    elseif D==3
        p=HEX
    else 
        @assert false
    end
    
    nodes, face_own_nodes = compute_nodes(p,orders)
    dofs = LagrangianDofBasis(T,nodes)
    reffaces = compute_lagrangian_reffaces(T,p,orders)
  
    nnodes = length(dofs.nodes)
    ndofs = length(dofs.dof_to_node)
    metadata = reffaces
    _reffaces = vcat(reffaces...)
    face_nodes = Gridap.ReferenceFEs._generate_face_nodes(nnodes,face_own_nodes,p,_reffaces)
    face_own_dofs = Gridap.ReferenceFEs._generate_face_own_dofs(face_own_nodes, dofs.node_and_comp_to_dof)
    face_dofs = Gridap.ReferenceFEs._generate_face_dofs(ndofs,face_own_dofs,p,_reffaces)
  
    if all(map(i->i==0,orders) ) && D>0
      conf = L2Conformity()
    else
      conf = GradConformity()
    end
  
    reffe = GenericRefFE{typeof(conf)}(
      ndofs,
      p,
      prebasis,
      dofs,
      conf,
      metadata,
      face_dofs)
  
    GenericLagrangianRefFE(reffe,face_nodes)
  end

  struct HQkIsoQ1Basis{T,D} <: AbstractVector{Field}
    order::Int
    function HQkIsoQ1Basis(::Type{T},D::Integer,order) where {T}
      @assert floor(log2(order)) == ceil(log2(order)) "The order of HQkIsoQ1Basis must be a power of 2"
      @assert D==1 || D==2
      new{T,D}(order)
    end
  end

  function Base.length(f::HQkIsoQ1Basis{T,D}) where {T,D}
    (f.order+1)^D
  end
  function Base.size(f::HQkIsoQ1Basis{T,D}) where {T,D}
    (length(f),)
  end
  # Dirty hack
  function Base.getindex(f::HQkIsoQ1Basis,i)
    # @warn "Calling a dirty hack implementation of HQkIsoQ1Basis"
    ConstantField(0.0)
  end
  
  # Methods required by RefFE re. the prebasis
  function Gridap.ReferenceFEs.get_orders(a::HQkIsoQ1Basis{T,D}) where {T,D}
    orders=zeros(Int,D)
    orders.=a.order
    Tuple(orders)
  end

  function Gridap.ReferenceFEs.get_order(a::HQkIsoQ1Basis{T,D}) where {T,D}
    a.order
  end


  function Gridap.Geometry.return_cache(a::HQkIsoQ1Basis{T,D},x::AbstractArray{<:Point}) where {T,D}
    num_levels = Int(floor(log2(a.order)))+1
    qk_iso_q1_reffes = Vector{GenericLagrangianRefFE}(undef, num_levels)
    for i=0:num_levels-1
      qk_iso_q1_reffes[i+1] = QkIsoQ1(T,D,2^i)
    end
    qk_iso_q1_reffes_at_x_cache = Vector{Any}(undef, num_levels)
    for i=1:num_levels
      shape_funs = get_shapefuns(qk_iso_q1_reffes[i])
      qk_iso_q1_reffes_at_x_cache[i] = return_cache(shape_funs, x)      
    end
    qk_iso_q1_reffes_at_x = Vector{Matrix{T}}(undef, num_levels)
    
    np   = length(x)
    ndof = num_dofs(qk_iso_q1_reffes[end]) 
    r = CachedArray(zeros(T,(np,ndof)))
    (r, qk_iso_q1_reffes,qk_iso_q1_reffes_at_x_cache,qk_iso_q1_reffes_at_x)
  end

  function Gridap.Geometry.evaluate!(cache,
                            a::HQkIsoQ1Basis{T,D},
                            x::AbstractArray{<:Point}) where {T,D}
    r, qk_iso_q1_reffes, qk_iso_q1_reffes_at_x_cache,qk_iso_q1_reffes_at_x = cache

    num_levels = length(qk_iso_q1_reffes)
    for i=1:num_levels
      shape_funs = get_shapefuns(qk_iso_q1_reffes[i])
      qk_iso_q1_reffes_at_x[i] = evaluate!(qk_iso_q1_reffes_at_x_cache[i],shape_funs, x)
    end

    np   = length(x)
    ndof = num_dofs(qk_iso_q1_reffes[end]) 
    setsize!(r, (np,ndof))

    order = a.order 
    
    # First vertices
    p = get_polytope(qk_iso_q1_reffes[1])
    for vertex=1:num_vertices(p)
      r[:,vertex]=qk_iso_q1_reffes_at_x[1][:,vertex]
    end

    for d=1:D 
      num_faces=length(get_faces(p,D,d)[1])
      for iface=1:num_faces
        last_lev_face_own_dofs_1=get_face_own_dofs(qk_iso_q1_reffes[end],1)[iface]
        last_lev_face_own_dofs_d=get_face_own_dofs(qk_iso_q1_reffes[end],d)[iface]
        for i=2:num_levels 
          k=2^(i-1)
          current_lev_face_own_dofs_1=get_face_own_dofs(qk_iso_q1_reffes[i],1)[iface]
          current_lev_face_dof = 1
          indices_current_lev = _generate_indices(current_lev_face_dof,
                                                  2,
                                                  length(current_lev_face_own_dofs_1))
          
          current_last_lev_face_dof = div(order,k)
          last_lev_face_dof_stride = current_last_lev_face_dof*2
          indices_last_lev    = _generate_indices(current_last_lev_face_dof,
                                                  last_lev_face_dof_stride,
                                                  length(last_lev_face_own_dofs_1))


          shape_current_lev=Vector{Int}(undef,d)
          shape_current_lev .= length(current_lev_face_own_dofs_1)
          shape_current_lev = Tuple(shape_current_lev)
          cis_current_lev=CartesianIndices(shape_current_lev)
          lis_current_lev=LinearIndices(shape_current_lev)

          shape_last_lev=Vector{Int}(undef,d)
          shape_last_lev .= length(last_lev_face_own_dofs_1)
          shape_last_lev = Tuple(shape_last_lev)
          cis_last_lev=CartesianIndices(shape_last_lev)
          lis_last_lev=LinearIndices(shape_last_lev)

          current_lev_face_own_dofs_d=get_face_own_dofs(qk_iso_q1_reffes[i],d)[iface]

          if d==1
            for (current_lev_face_dof,current_last_lev_face_dof) in zip(indices_current_lev,indices_last_lev)
              @debug "$(d) $(i) $(iface): $(current_lev_face_dof) $(current_last_lev_face_dof)"
              r[:,last_lev_face_own_dofs_1[current_last_lev_face_dof]]=
                    qk_iso_q1_reffes_at_x[i][:,current_lev_face_own_dofs_1[current_lev_face_dof]]
             end
          elseif d==2
            indices_current_lev_full = _generate_indices(current_lev_face_dof,
                                                         1,
                                                         length(current_lev_face_own_dofs_1))
                                                      
            indices_last_lev_full = _generate_indices(current_last_lev_face_dof,
                                                      div(last_lev_face_dof_stride,2),
                                                      length(last_lev_face_own_dofs_1))
                                        
            for (i1,i2) in zip(indices_current_lev_full,indices_last_lev_full)
              if mod(i1,2)!=0 
                a = indices_current_lev_full
                b = indices_last_lev_full
              else
                a = indices_current_lev
                b = indices_last_lev
              end
              for (j1,j2) in zip(a,b)
                li1=lis_current_lev[CartesianIndex(i1,j1)]
                li2=lis_last_lev[CartesianIndex(i2,j2)]
                @debug "$(d) $(i) $(iface): $(li1) $(li2)"
                r[:,last_lev_face_own_dofs_d[li2]]=
                   qk_iso_q1_reffes_at_x[i][:,current_lev_face_own_dofs_d[li1]]
              end 
            end
          end 
        end
      end
    end
    r.array
  end 

  function _generate_indices(start,stride,n)
    collect(start:stride:n)
  end 
  
  struct LinearCombinationDofBasis{T,B<:AbstractMatrix} <: AbstractVector{T}
    basis::AbstractVector{T}
    matrix::B
    function LinearCombinationDofBasis(basis::A,matrix::B) where {A,B}
      @assert length(basis)==size(matrix,2)
      T = eltype(basis)
      @assert T <: Dof
      new{T,B}(basis,matrix)
    end
  end 

  function Base.length(b::LinearCombinationDofBasis)
    size(b.matrix,1)
  end

  function Base.size(b::LinearCombinationDofBasis)
    (length(b),)
  end

  function Base.getindex(b::LinearCombinationDofBasis,i)
    @notimplemented 
  end

  # We need to implement this to be able to use the basis with the skin of a Lagrangian basis
  function Base.getproperty(b::LinearCombinationDofBasis,i::Symbol)
    if i==:node_and_comp_to_dof
      b.basis.node_and_comp_to_dof
    elseif i==:dof_to_node
      b.basis.dof_to_node
    elseif i==:dof_to_comp
      b.basis.dof_to_comp
    elseif i==:nodes
      b.basis.nodes
    else
      getfield(b,i)
    end
  end

  function Gridap.Geometry.return_cache(b::LinearCombinationDofBasis,field)
    cbasis = return_cache(b.basis,field)
    T=eltype(b.matrix)
    v = copy(evaluate!(cbasis,b.basis,field))
    (cbasis,v)
  end

  using LinearAlgebra
  function Gridap.Geometry.evaluate!(cache,b::LinearCombinationDofBasis,field)
    cbasis,v = cache
    b_basis=evaluate!(cbasis,b.basis,field)
    LinearAlgebra.mul!(v,b.matrix,b_basis)
    v
  end

  function HQkIsoQ1(::Type{T},D::Integer,order) where {T}
    @assert D==1 || D==2
    @assert floor(log2(order)) == ceil(log2(order)) "The order of the Lagrangian reference FE must be a power of 2"
    qk_iso_q1_reffe = QkIsoQ1(T,D,order)
    hierarchical_basis=HQkIsoQ1Basis(T,D,order)
    x=get_nodes(qk_iso_q1_reffe.reffe.dofs)
    W=evaluate(hierarchical_basis,x)
    hierarchical_dof_basis=LinearCombinationDofBasis(qk_iso_q1_reffe.reffe.dofs,inv(W))
    reffe = GenericRefFE{typeof(qk_iso_q1_reffe.reffe.conformity)}(
      qk_iso_q1_reffe.reffe.ndofs,
      qk_iso_q1_reffe.reffe.polytope,
      hierarchical_basis,
      hierarchical_dof_basis,
      qk_iso_q1_reffe.reffe.conformity,
      qk_iso_q1_reffe.reffe.metadata,
      qk_iso_q1_reffe.reffe.face_dofs)

    GenericLagrangianRefFE(reffe,qk_iso_q1_reffe.reffe.face_dofs)
  end

  # function _compute_high_order_nodes(p::Polytope{D},orders) where D
  #   nodes = Point{D,Float64}[]
  #   facenodes = [Int[] for i in 1:num_faces(p)]
  #   _compute_high_order_nodes_dim_0!(nodes,facenodes,p)
  #   for d in 1:(num_dims(p)-1)
  #     _compute_high_order_nodes_dim_d!(nodes,facenodes,p,orders,Val{d}())
  #   end
  #   _compute_high_order_nodes_dim_D!(nodes,facenodes,p,orders)
  #   (nodes, facenodes)
  # end
  
  # function _compute_high_order_nodes_dim_0!(nodes,facenodes,p)
  #   x = get_vertex_coordinates(p)
  #   k = 1
  #   for vertex in 1:num_vertices(p)
  #     push!(nodes,x[vertex])
  #     push!(facenodes[vertex],k)
  #     k += 1
  #   end
  # end
  
  # @noinline function _compute_high_order_nodes_dim_d!(nodes,facenodes,p,orders,::Val{d}) where d
  #   x = get_vertex_coordinates(p)
  #   offset = get_offset(p,d)
  #   k = length(nodes)+1
  #   for iface in 1:num_faces(p,d)
  #     face = Polytope{d}(p,iface)
  #     face_ref_x = get_vertex_coordinates(face)
  #     face_prebasis = MonomialBasis(Float64,face,1)
  #     change = inv(evaluate(face_prebasis,face_ref_x))
  #     face_shapefuns = linear_combination(change,face_prebasis)
  #     face_vertex_ids = get_faces(p,d,0)[iface]
  #     face_x = x[face_vertex_ids]
  #     face_orders = compute_face_orders(p,face,iface,orders)
  #     face_interior_nodes = compute_own_nodes(face,face_orders)
  #     face_high_x = evaluate(face_shapefuns,face_interior_nodes)*face_x
  #     for xi in 1:length(face_high_x)
  #       push!(nodes,face_high_x[xi])
  #       push!(facenodes[iface+offset],k)
  #       k += 1
  #     end
  #   end
  # end
  
  # function _compute_high_order_nodes_dim_D!(nodes,facenodes,p,orders)
  #   k = length(nodes)+1
  #   p_high_x = compute_own_nodes(p,orders)
  #   for xi in 1:length(p_high_x)
  #     push!(nodes,p_high_x[xi])
  #     push!(facenodes[end],k)
  #     k += 1
  #   end
  # end