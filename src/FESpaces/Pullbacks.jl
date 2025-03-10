
# TODO: We probably want to export these from Gridap.ReferenceFEs
using Gridap.ReferenceFEs: PushforwardRefFE, Pushforward, Pullback
using Gridap.ReferenceFEs: ContraVariantPiolaMap, CoVariantPiolaMap

function get_cell_dof_basis(
  model::DiscreteModel, cell_reffe::AbstractArray{<:GenericRefFE{T}}, conformity::Conformity
) where T <: PushforwardRefFE
  pushforward, cell_change, cell_args = get_cell_pushforward(
    Pushforward(T), model, cell_reffe, conformity
  )
  cell_ref_dofs = lazy_map(get_dof_basis, cell_reffe)
  cell_phy_dofs = lazy_map(inverse_map(Pullback(pushforward)), cell_ref_dofs, cell_args...)
  return lazy_map(linear_combination, cell_change, cell_phy_dofs) # TODO: Inverse and transpose
end

function get_cell_shapefuns(
  model::DiscreteModel, cell_reffe::AbstractArray{<:GenericRefFE{T}}, conformity::Conformity
) where T <: PushforwardRefFE
  pushforward, cell_change, cell_args = get_cell_pushforward(
    Pushforward(T), model, cell_reffe, conformity
  )
  cell_ref_fields = lazy_map(get_shapefuns, cell_reffe)
  cell_phy_fields = lazy_map(pushforward, cell_ref_fields, cell_args...)
  return lazy_map(linear_combination, cell_change, cell_phy_fields)
end

function get_cell_pushforward(
  ::Pushforward, model::DiscreteModel, cell_reffe, conformity
)
  @abstractmethod
end

# ContraVariantPiolaMap

function get_cell_pushforward(
  ::ContraVariantPiolaMap, model::DiscreteModel, cell_reffe, conformity
)
  cell_map  = get_cell_map(get_grid(model))
  Jt = lazy_map(Broadcasting(∇),cell_map)
  change = get_sign_flip(model, cell_reffe)
  return ContraVariantPiolaMap(), change, (Jt,)
end

# CoVariantPiolaMap

function get_cell_pushforward(
  ::CoVariantPiolaMap, model::DiscreteModel, cell_reffe, conformity
)
  cell_map = get_cell_map(get_grid(model))
  Jt = lazy_map(Broadcasting(∇),cell_map)
  change = lazy_map(r -> Diagonal(ones(num_dofs(r))), cell_reffe) # TODO: Replace by edge-signs
  return CoVariantPiolaMap(), change, (Jt,)
end

# NormalSignMap

"""
    struct NormalSignMap <: Map
      ...
    end
"""
struct NormalSignMap{T} <: Map
  model::T
  facet_owners::Vector{Int32}
end

function NormalSignMap(model)
  facet_owners = compute_facet_owners(model)
  NormalSignMap(model,facet_owners)
end

function return_value(k::NormalSignMap,reffe,facet_own_dofs,cell)
  Diagonal(fill(one(Float64), num_dofs(reffe)))
end

function return_cache(k::NormalSignMap,reffe,facet_own_dofs,cell)
  model = k.model
  Dc = num_cell_dims(model)
  topo = get_grid_topology(model)

  cell_facets = get_faces(topo, Dc, Dc-1)
  cell_facets_cache = array_cache(cell_facets)

  return cell_facets, cell_facets_cache, CachedVector(Float64)
end

function evaluate!(cache,k::NormalSignMap,reffe,facet_own_dofs,cell)
  cell_facets,cell_facets_cache,dof_sign_cache = cache
  facet_owners = k.facet_owners

  setsize!(dof_sign_cache, (num_dofs(reffe),))
  dof_sign = dof_sign_cache.array
  fill!(dof_sign, one(eltype(dof_sign)))

  facets = getindex!(cell_facets_cache,cell_facets,cell)
  for (lfacet,facet) in enumerate(facets)
    owner = facet_owners[facet]
    if owner != cell
      for dof in facet_own_dofs[lfacet]
        dof_sign[dof] = -1.0
      end
    end
  end

  return Diagonal(dof_sign)
end

function get_sign_flip(model::DiscreteModel{Dc}, cell_reffe) where Dc
  # Comment: lazy_maps on cell_reffes are very optimised, since they are CompressedArray/FillArray
  get_facet_own_dofs(reffe) = view(get_face_own_dofs(reffe),get_dimrange(get_polytope(reffe),Dc-1))
  cell_facet_own_dofs = lazy_map(get_facet_own_dofs, cell_reffe)
  cell_ids = IdentityVector(Int32(num_cells(model)))
  return lazy_map(NormalSignMap(model), cell_reffe, cell_facet_own_dofs, cell_ids)
end

function compute_facet_owners(model::DiscreteModel{Dc}) where {Dc}
  topo = get_grid_topology(model)
  facet_to_cell = get_faces(topo, Dc-1, Dc)

  nfacets = num_faces(topo, Dc-1)
  owners = Vector{Int32}(undef, nfacets)
  for facet in 1:nfacets
    facet_cells = view(facet_to_cell, facet)
    owners[facet] = first(facet_cells)
  end

  return owners
end
