
# TODO: We probably want to export these from Gridap.ReferenceFEs
using Gridap.ReferenceFEs: PushforwardRefFE, Pushforward, InversePullback
using Gridap.ReferenceFEs: ContraVariantPiolaMap, CoVariantPiolaMap

function get_cell_dof_basis(
  model::DiscreteModel, cell_reffe::AbstractArray{<:GenericRefFE{T}}, conformity::Conformity
) where T <: PushforwardRefFE
  pushforward = Pushforward(T)
  cell_args = get_cell_pushforward_arguments(pushforward, model, cell_reffe, conformity)
  cell_dofs = lazy_map(get_dof_basis, cell_reffe)
  return lazy_map(InversePullback(pushforward), cell_dofs, cell_args...)
end

function get_cell_shapefuns(
  model::DiscreteModel, cell_reffe::AbstractArray{<:GenericRefFE{T}}, conformity::Conformity
) where T <: PushforwardRefFE
  pushforward = Pushforward(T)
  cell_args = get_cell_pushforward_arguments(pushforward, model, cell_reffe, conformity)
  cell_dofs = lazy_map(get_shapefuns, cell_reffe)
  return lazy_map(pushforward, cell_dofs, cell_args...)
end

function get_cell_pushforward_arguments(
  ::Pushforward, model::DiscreteModel, cell_reffe, conformity
)
  @abstractmethod
end

# ContraVariantPiolaMap

function get_cell_pushforward_arguments(
  ::ContraVariantPiolaMap, model::DiscreteModel, cell_reffe, conformity
)
  cell_map  = get_cell_map(get_grid(model))
  Jt = lazy_map(Broadcasting(∇),cell_map)
  sign_flip = lazy_map(Broadcasting(constant_field),get_sign_flip(model, cell_reffe))
  return Jt, sign_flip
end

# CoVariantPiolaMap

function get_cell_pushforward_arguments(
  ::CoVariantPiolaMap, model::DiscreteModel, cell_reffe, conformity
)
  cell_map  = get_cell_map(get_grid(model))
  Jt = lazy_map(Broadcasting(∇),cell_map)
  return Jt
end

# SignFlipMap

struct SignFlipMap{T} <: Map
  model::T
  facet_owners::Vector{Int32}
end

function SignFlipMap(model)
  facet_owners = compute_facet_owners(model)
  SignFlipMap(model,facet_owners)
end

function return_value(k::SignFlipMap,reffe,facet_own_dofs,cell)
  fill(false, num_dofs(reffe))
end

function return_cache(k::SignFlipMap,reffe,facet_own_dofs,cell)
  model = k.model
  Dc = num_cell_dims(model)
  topo = get_grid_topology(model)

  cell_facets = get_faces(topo, Dc, Dc-1)
  cell_facets_cache = array_cache(cell_facets)

  return cell_facets, cell_facets_cache, CachedVector(Bool)
end

function evaluate!(cache,k::SignFlipMap,reffe,facet_own_dofs,cell)
  cell_facets,cell_facets_cache,sign_flip_cache = cache
  facet_owners = k.facet_owners

  setsize!(sign_flip_cache, (num_dofs(reffe),))
  sign_flip = sign_flip_cache.array
  sign_flip .= false

  facets = getindex!(cell_facets_cache,cell_facets,cell)
  for (lfacet,facet) in enumerate(facets)
    owner = facet_owners[facet]
    if owner != cell
      for dof in facet_own_dofs[lfacet]
        sign_flip[dof] = true
      end
    end
  end

  return sign_flip
end

function get_sign_flip(model::DiscreteModel{Dc}, cell_reffe) where Dc
  # Comment: lazy_maps on cell_reffes are very optimised, since they are CompressedArray/FillArray
  get_facet_own_dofs(reffe) = view(get_face_own_dofs(reffe),get_dimrange(get_polytope(reffe),Dc-1))
  cell_facet_own_dofs = lazy_map(get_facet_own_dofs, cell_reffe)
  cell_ids = IdentityVector(Int32(num_cells(model)))
  return lazy_map(SignFlipMap(model), cell_reffe, cell_facet_own_dofs, cell_ids)
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
