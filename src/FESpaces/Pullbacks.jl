
function get_cell_dof_basis(
  model::DiscreteModel, cell_reffe::AbstractArray{T}, conf::Conformity
) where T <: ReferenceFE

  pushforward = Pushforward(get_name(T), conf)
  if pushforward isa IdentityPiolaMap
    lazy_map(get_dof_basis,cell_reffe)
  else
    pushforward, cell_change, cell_args = get_cell_pushforward(
      pushforward, model, cell_reffe)
    cell_ref_dofs = lazy_map(get_dof_basis, cell_reffe)
    cell_phy_dofs = lazy_map(inverse_map(Pullback(pushforward)), cell_ref_dofs, cell_args...)
    lazy_map(linear_combination, cell_change, cell_phy_dofs) # TODO: Inverse and transpose
  end
end

function get_cell_shapefuns(
  model::DiscreteModel, cell_reffe::AbstractArray{T}, conf::Conformity
) where T <: ReferenceFE

  pushforward = Pushforward(get_name(T), conf)
  if pushforward isa IdentityPiolaMap
    lazy_map(get_shapefuns,cell_reffe)
  else
    pushforward, cell_change, cell_args = get_cell_pushforward(
      pushforward, model, cell_reffe)
    cell_ref_fields = lazy_map(get_shapefuns, cell_reffe)
    cell_phy_fields = lazy_map(pushforward, cell_ref_fields, cell_args...)
    lazy_map(linear_combination, cell_change, cell_phy_fields)
  end
end

function get_cell_pushforward(
  ::Pushforward, model::DiscreteModel, cell_reffe
)
  @abstractmethod
end

# ContraVariantPiolaMap

function get_cell_pushforward(
  p::ContraVariantPiolaMap, model::DiscreteModel, cell_reffe
)
  cell_map  = get_cell_map(get_grid(model))
  Jt = lazy_map(Broadcasting(∇),cell_map)
  change = get_sign_flip(model, cell_reffe)
  return p, change, (Jt,)
end

# CoVariantPiolaMap

function get_cell_pushforward(
  p::CoVariantPiolaMap, model::DiscreteModel, cell_reffe
)
  cell_map = get_cell_map(get_grid(model))
  Jt = lazy_map(Broadcasting(∇),cell_map)
  change = lazy_map(r -> Diagonal(ones(num_dofs(r))), cell_reffe) # TODO: Replace by edge-signs
  return p, change, (Jt,)
end

# DoubleContraVariantPiolaMap

using Gridap.ReferenceFEs: DoubleContraVariantPiolaMap
function get_cell_pushforward(
  ::DoubleContraVariantPiolaMap, model::DiscreteModel, cell_reffe, conformity
)
  cell_map = get_cell_map(get_grid(model))
  Jt = lazy_map(Broadcasting(∇),cell_map)
  change = lazy_map(r -> Diagonal(fill(one(Float64), num_dofs(r))), cell_reffe)
  #change = get_sign_flip(model, cell_reffe)
  return DoubleContraVariantPiolaMap(), change, (Jt,)
end

# NormalSignMap

"""
    struct NormalSignMap <: Map
      ...
    end

The `NormalSignMap` compute the signs to apply to the mapped reference normals,
for each facet of a physical cell.

Each physical facet ``f`` is shared by up to two cells ``K`` and ``K'``. The
orientation of the physical/global normal to ``f`` is chosen by the main cell
``K``, the first one in the list of adjascent cells to ``f`` in the grid topology.

The physical/global normal is the (normalized) Piola mapped reference normal to
``f̂ = F⁻¹(f)`` where ``F`` is the geometrical map ``F:K̂->K``. It is also minus
the (normalized) Piola mapped reference normal to ``f̂' = F'⁻¹(f)`` where ``F'``
is the geometrical map ``F':K̂->K'``.
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

  o = one(eltype(dof_sign))
  fill!(dof_sign, o)

  facets = getindex!(cell_facets_cache,cell_facets,cell)
  for (lfacet,facet) in enumerate(facets)
    owner = facet_owners[facet]
    if owner != cell
      for dof in facet_own_dofs[lfacet]
        dof_sign[dof] = -o
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
