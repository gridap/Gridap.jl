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