
struct PolytopalFESpace{V,M} <: SingleFieldFESpace
  vector_type::Type{V}
  ndofs::Int
  cell_dofs_ids::AbstractArray
  fe_basis::CellField
  metadata::M
end

# Constructors 

function PolytopalFESpace(
  trian::Triangulation{Dc}, order::Integer; T = Float64, vector_type = Vector{T}
) where Dc
  ncells = num_cells(trian)
  prebasis = Polynomials.MonomialBasis{Dc}(T, order)
  cell_prebasis = Fill(prebasis, ncells)
  cell_change = lazy_map(centroid_map, get_polytopes(trian))
  cell_shapefuns = lazy_map(Broadcasting(∘), cell_prebasis, cell_change)
  fe_basis = SingleFieldFEBasis(cell_shapefuns, trian, TestBasis(), PhysicalDomain())
  cell_dof_ids, ndofs = compute_discontinuous_cell_dofs(
    Base.OneTo(ncells), Fill(length(prebasis), ncells)
  )
  return PolytopalFESpace(
    vector_type,ndofs,cell_dof_ids,fe_basis,nothing
  )
end

function PolytopalFESpace(model::DiscreteModel, args...; kwargs...)
  PolytopalFESpace(Triangulation(model), args...; kwargs...)
end

# FESpace interface

ConstraintStyle(::Type{<:PolytopalFESpace}) = UnConstrained()
get_free_dof_ids(f::PolytopalFESpace) = Base.OneTo(f.ndofs)
get_fe_basis(f::PolytopalFESpace) = f.fe_basis
get_fe_dof_basis(f::PolytopalFESpace) = @notimplemented
get_cell_dof_ids(f::PolytopalFESpace) = f.cell_dofs_ids
get_triangulation(f::PolytopalFESpace) = get_triangulation(f.fe_basis)
get_dof_value_type(f::PolytopalFESpace{V}) where V = eltype(V)
get_vector_type(f::PolytopalFESpace{V}) where V = V
get_cell_is_dirichlet(f::PolytopalFESpace) = Fill(false,num_cells(get_triangulation(f)))

# SingleFieldFESpace interface

get_dirichlet_dof_ids(f::PolytopalFESpace) = Base.OneTo(0)
num_dirichlet_tags(f::PolytopalFESpace) = 0
get_dirichlet_dof_tag(f::PolytopalFESpace) = Int8[]

function scatter_free_and_dirichlet_values(f::PolytopalFESpace,free_values,dirichlet_values)
  @check length(dirichlet_values) == 0 "PolytopalFESpace does not support Dirichlet values"
  cell_dof_ids = get_cell_dof_ids(f)
  lazy_map(Broadcasting(Reindex(free_values)),cell_dof_ids)
end

function gather_free_and_dirichlet_values!(free_vals,dirichlet_vals,f::PolytopalFESpace,cell_vals)
  cell_dofs = get_cell_dof_ids(f)
  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)
  cells = 1:length(cell_vals)
  _free_and_dirichlet_values_fill!(
    free_vals, dirichlet_vals,
    cache_vals, cache_dofs,
    cell_vals, cell_dofs, cells
  )
  return free_vals, dirichlet_vals
end

function gather_dirichlet_values!(dirichlet_vals,f::PolytopalFESpace,cell_vals)
  @check length(dirichlet_vals) == 0 "PolytopalFESpace does not support Dirichlet values"
  dirichlet_vals
end

# Change of coordinate map

struct CentroidCoordinateChangeMap <: Map end

function centroid_map(poly::Polytope{D}) where D
  v = get_vertex_coordinates(poly)
  xc = mean(v)
  h = maximum(vk -> norm(vk-xc), v)
  origin = xc/h
  gradient = TensorValues.diagonal_tensor(VectorValue([h for i in 1:D]))
  return affine_map(gradient,origin)
end

function Arrays.lazy_map(::typeof(evaluate),a::LazyArray{<:Fill{typeof(centroid_map)}},x::AbstractVector{<:Point})
  polys = a.args[1]
  lazy_map(CentroidCoordinateChangeMap(),polys,x)
end

function Arrays.evaluate!(cache,::CentroidCoordinateChangeMap,poly::Polytope{D},x::Point{D}) where D
  v = get_vertex_coordinates(poly)
  xc = mean(v)
  h = maximum(vk -> norm(vk-xc), v)
  return (x - xc) / h
end

function Arrays.return_cache(::CentroidCoordinateChangeMap,poly::Polytope{D},x::AbstractVector{Point{D}}) where D
  return CachedArray(similar(x))
end

function Arrays.evaluate!(cache,::CentroidCoordinateChangeMap,poly::Polytope{D},x::AbstractVector{Point{D}}) where D
  setsize!(cache,length(x))
  y = cache.array
  v = get_vertex_coordinates(poly)
  xc = mean(v)
  h = maximum(vk -> norm(vk-xc), v)
  for i in eachindex(x)
    y[i] = (x[i] - xc) / h
  end
  return y
end
