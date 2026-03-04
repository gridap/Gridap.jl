
struct PolytopalFESpace{V,M} <: SingleFieldFESpace
  vector_type::Type{V}
  nfree::Int
  ndirichlet::Int
  cell_dofs_ids::AbstractArray
  fe_basis::CellField
  cell_is_dirichlet::Vector{Bool}
  dirichlet_dof_tag::Vector{Int8}
  dirichlet_cells::Vector{Int32}
  ntags::Int
  order::Int
  metadata::M
end

# Constructors 

function _filter_from_space(space::Symbol)
  if space == :P
    return Polynomials._p_filter
  elseif space == :Q
    return Polynomials._q_filter
  elseif space == :S
    return Polynomials._s_filter
  else
    @notimplemented
  end
end

function _prebasis_from_space(D,T,order,space)
  if space == :RT
    @assert D >= 2
    prebasis = Polynomials.PCurlGradMonomialBasis{D}(T, order)
  elseif space == :ND
    @assert D >= 2
    prebasis = Polynomials.NedelecPrebasisOnSimplex{D}(order)
  else
    prebasis = Polynomials.MonomialBasis{D}(T, order, _filter_from_space(space))
  end
  return prebasis
end

function PolytopalFESpace(
  t::Triangulation,args...;trian=nothing,kwargs...
)
  @check isnothing(trian)
  model = get_active_model(t)
  PolytopalFESpace(model,args...;trian=t,kwargs...)
end

function PolytopalFESpace(
  model::DiscreteModel,::Type{T},order::Integer;
  space=:P,vector_type=nothing,hierarchical=false,
  kwargs...
) where {T}
  D = num_cell_dims(model)
  prebasis = _prebasis_from_space(D,T,order,space)
  if hierarchical
    lt(t1,t2) = (maximum(Tuple(t1)) < maximum(Tuple(t2))) || isless(t1,t2)
    sort!(prebasis.terms;lt)
  end
  cell_prebasis = Fill(prebasis, num_cells(model))
  vtype = ifelse(!isnothing(vector_type),vector_type,Vector{_dof_type(T)})
  PolytopalFESpace(vtype,model,cell_prebasis;kwargs...)
end

function PolytopalFESpace(
  vector_type::Type,
  model::DiscreteModel,
  cell_prebasis::AbstractArray; 
  trian = Triangulation(model),
  labels = get_face_labeling(model),
  dirichlet_tags = Int[],
  dirichlet_masks = nothing,
  orthonormal = false,
  local_kernel = nothing
)
  T = return_type(first(cell_prebasis))
  Dc = num_cell_dims(model)

  order = maximum(Polynomials.get_order,cell_prebasis)
  cell_polytopes = Geometry.get_cell_polytopes(model)
  if eltype(cell_polytopes) <: ReferenceFEs.GeneralPolytope
    cell_change = lazy_map(centroid_map, cell_polytopes)
    cell_shapefuns = lazy_map(Broadcasting(∘), cell_prebasis, cell_change)
    domain_style = PhysicalDomain()
  else
    cell_shapefuns = cell_prebasis
    domain_style = ReferenceDomain()
  end
  if !isnothing(local_kernel)
    cell_shapefuns = remove_local_kernel(cell_shapefuns, trian, T, order, local_kernel)
  end
  if orthonormal
    cell_shapefuns = orthogonalise_basis(cell_shapefuns, trian, order)
  end
  fe_basis = SingleFieldFEBasis(cell_shapefuns, trian, TestBasis(), domain_style)
  cell_conformity = DiscontinuousCellConformity(cell_shapefuns)

  ntags = length(dirichlet_tags)
  if ntags != 0
    @notimplementedif !isnothing(local_kernel)
    cell_to_tag = get_face_tag_index(labels,dirichlet_tags,Dc)
    cell_is_dirichlet = map(!isequal(UNSET),cell_to_tag)
    cell_dof_ids, nfree, ndir, dirichlet_dof_tag, dirichlet_cells = compute_discontinuous_cell_dofs(
      cell_conformity, cell_to_tag, dirichlet_masks
    )
  else
    ndir = 0
    dirichlet_dof_tag = Int8[]
    dirichlet_cells = Int32[]
    cell_is_dirichlet = fill(false,num_cells(trian))
    cell_dof_ids, nfree = compute_discontinuous_cell_dofs(cell_conformity)
  end

  metadata = cell_conformity
  return PolytopalFESpace(
    vector_type,nfree,ndir,cell_dof_ids,fe_basis,
    cell_is_dirichlet,dirichlet_dof_tag,dirichlet_cells,ntags,
    order,metadata
  )
end

# FESpace interface

ConstraintStyle(::Type{<:PolytopalFESpace}) = UnConstrained()
get_free_dof_ids(f::PolytopalFESpace) = Base.OneTo(f.nfree)
get_fe_basis(f::PolytopalFESpace) = f.fe_basis
get_fe_dof_basis(f::PolytopalFESpace) = @notimplemented
get_cell_dof_ids(f::PolytopalFESpace) = f.cell_dofs_ids
get_triangulation(f::PolytopalFESpace) = get_triangulation(f.fe_basis)
get_dof_value_type(f::PolytopalFESpace{V}) where V = eltype(V)
get_vector_type(f::PolytopalFESpace{V}) where V = V
get_cell_is_dirichlet(f::PolytopalFESpace) = f.cell_is_dirichlet
get_cell_conformity(f::PolytopalFESpace{V,<:CellConformity}) where V = f.metadata

# SingleFieldFESpace interface

get_dirichlet_dof_ids(f::PolytopalFESpace) = Base.OneTo(f.ndirichlet)
num_dirichlet_tags(f::PolytopalFESpace) = f.ntags
get_dirichlet_dof_tag(f::PolytopalFESpace) = f.dirichlet_dof_tag

function _cell_vals(f::PolytopalFESpace,object)
  dΩ = Measure(get_triangulation(f),2*f.order)
  u, v = get_trial_fe_basis(f), get_fe_basis(f)

  lhs = get_array(∫(u⋅v)dΩ)
  rhs = get_array(∫(v⋅object)dΩ)
  cell_vals = lazy_map(LocalSolveMap(),lhs,rhs)
  return cell_vals
end

_cell_vals(f::TrialFESpace{<:PolytopalFESpace},object) = _cell_vals(f.space,object)

function scatter_free_and_dirichlet_values(f::PolytopalFESpace,free_values,dirichlet_values)
  @check eltype(free_values) == eltype(dirichlet_values) """\n
  The entries stored in free_values and dirichlet_values should be of the same type.

  This error shows up e.g. when trying to build a FEFunction from a vector of integers
  if the Dirichlet values of the underlying space are of type Float64, or when the
  given free values are Float64 and the Dirichlet values ComplexF64.
  """
  cell_dof_ids = get_cell_dof_ids(f)
  lazy_map(Broadcasting(PosNegReindex(free_values,dirichlet_values)),cell_dof_ids)
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
  cell_dofs = get_cell_dof_ids(f)
  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)
  free_vals = zero_free_values(f)
  cells = f.dirichlet_cells
  _free_and_dirichlet_values_fill!(
    free_vals, dirichlet_vals,
    cache_vals, cache_dofs,
    cell_vals, cell_dofs, cells
  )
  return dirichlet_vals
end

# Change of coordinate map

function centroid_map(poly::Polytope)
  D, Dp = num_dims(poly), num_point_dims(poly)
  centroid_map(Val(D),Val(Dp),get_vertex_coordinates(poly))
end

function centroid_map(::Val{D},::Val{D},vertices) where D
  pmin, pmax = get_bounding_box(vertices)
  xc = 0.5 * (pmin + pmax)
  h = 0.5 * (pmax - pmin)
  G = TensorValues.diagonal_tensor(VectorValue(ntuple(i -> 1.0/h[i], Val(D))))
  O = -xc⋅G
  return affine_map(G,O)
end

function centroid_map(::Val{D},::Val{Dp},vertices) where {D,Dp}
  @check Dp > D
  t_space = ReferenceFEs.compute_tangent_space(Val(D),vertices)
  K = TensorValues.tensor_from_columns(t_space...)
  pmin, pmax = get_bounding_box(v -> v⋅K, vertices)
  xc = 0.5 * (pmin + pmax)
  h = 0.5 * (pmax - pmin)
  G = TensorValues.diagonal_tensor(VectorValue(ntuple(i -> 1.0/h[i], Val(D))))
  O = -xc⋅G
  return affine_map(K⋅G,O)
end

# Orthogonalisation

function orthogonalise_basis(cell_basis,trian,order)
  cell_quads = Quadrature(trian,2*order)
  cell_integ = lazy_map(Broadcasting(Operation(⋅)),cell_basis,lazy_map(transpose,cell_basis))
  cell_vals  = lazy_map(evaluate,cell_integ,lazy_map(get_coordinates,cell_quads))
  cell_mass  = lazy_map(IntegrationMap(),cell_vals,lazy_map(get_weights,cell_quads))
  cell_coeffs = lazy_map(OrthogonaliseBasisMap(),cell_mass)
  return lazy_map(linear_combination,cell_coeffs,cell_basis)
end

struct OrthogonaliseBasisMap <: Map end

function Arrays.return_cache(::OrthogonaliseBasisMap,M::Matrix)
  return CachedArray(similar(M))
end

function Arrays.evaluate!(cache,::OrthogonaliseBasisMap,M::Matrix)
  setsize!(cache,size(M))
  N = cache.array
  gram_shmidt!(N,M)
  return N
end

# Orthogonalise a basis against itself, with respect 
# of the inner product defined by M
# An alternative would be N = inv(cholesky(M).U)
function gram_shmidt!(N,M)
  n = size(M,1)
  fill!(N,0.0)
  for i in 1:n
    N[i,i] = 1.0
  end

  Nk = eachcol(N)
  for i in 1:n
    Ni = Nk[i]
    for j in 1:i-1
      Nj = Nk[j]
      rij = dot(Ni, M, Nj)
      Ni .-= rij .* Nj
    end
    rii = sqrt(dot(Ni, M , Ni))
    Ni ./= rii
  end

  return N
end

# Local kernel removal 

function remove_local_kernel(cell_basis,trian,T,order,local_kernel)
  if isa(local_kernel,Function) 
    local_kernel_func = local_kernel
  else 
    local_kernel_func = _kernel_from_symbol(local_kernel,T,cell_basis)
  end
  
  cell_quads = Quadrature(trian,2*order)
  cell_integ = local_kernel_func(lazy_map(transpose,cell_basis))
  cell_vals = lazy_map(evaluate,cell_integ,lazy_map(get_coordinates,cell_quads))
  cell_kernel = lazy_map(IntegrationMap(),cell_vals,lazy_map(get_weights,cell_quads))
  cell_coeffs = lazy_map(NullspaceMap(),cell_kernel)
  return lazy_map(linear_combination,cell_coeffs,cell_basis)
end

# Stolen from the MomentBased branch
component_basis(T::Type{<:Real}) = [one(T)]
function component_basis(V::Type{<:MultiValue})
  T = eltype(V)
  n = num_components(V)
  z, o = zero(T), one(T)
  return [V(ntuple(i -> ifelse(i == j, o, z),Val(n))) for j in 1:n]
end

function _kernel_from_symbol(k::Symbol,T,cell_basis)
  if k == :constants
    fields = map(constant_field,component_basis(T))
    cell_fields = Fill(fields, length(cell_basis))
    kernel(cell_basis) = lazy_map(Broadcasting(Operation(⊙)),cell_fields,cell_basis)
  else
    @notimplemented
  end
  return kernel
end

struct NullspaceMap{T} <: Map
  tol::T
end

NullspaceMap() = NullspaceMap(10*eps(Float64))

function Arrays.return_cache(::NullspaceMap,K::Matrix)
  return CachedArray(similar(K))
end

function Arrays.evaluate!(cache,k::NullspaceMap,K::Matrix)
  f = svd(K;full=true) # If K is not square, there svd! doesn't really use cache
  
  n = size(K, 2)
  m = sum(s -> s > k.tol, f.S)
  setsize!(cache, (n, n - m))
  N = cache.array

  for i in 1:n
    for j in (m+1):n
      N[i,j-m] = f.Vt[j,i]
    end
  end

  return N
end

# struct CentroidCoordinateChangeMap <: Map end
# 
# function Arrays.lazy_map(::typeof(evaluate),a::LazyArray{<:Fill{typeof(centroid_map)}},x::AbstractVector)
#   polys = a.args[1]
#   lazy_map(CentroidCoordinateChangeMap(),polys,x)
# end
# 
# function Arrays.evaluate!(cache,::CentroidCoordinateChangeMap,poly::Polytope,x::Point)
#   pmin, pmax = get_bounding_box(poly)
#   xc = 0.5 * (pmin + pmax)
#   h = 0.5 * (pmax - pmin)
#   return (x - xc) ./ h
# end
# 
# function Arrays.return_cache(::CentroidCoordinateChangeMap,poly::Polytope,x::AbstractVector{<:Point})
#   return CachedArray(similar(x))
# end
# 
# function Arrays.evaluate!(cache,::CentroidCoordinateChangeMap,poly::Polytope,x::AbstractVector{<:Point})
#   setsize!(cache,size(x))
#   y = cache.array
#   pmin, pmax = get_bounding_box(poly)
#   xc = 0.5 * (pmin + pmax)
#   h = 0.5 * (pmax - pmin)
#   for i in eachindex(x)
#     y[i] = (x[i] - xc) ./ h
#   end
#   return y
# end

##################

function shoelace(face_ents)
  shift = circshift(face_ents, -1)
  area_components = map(face_ents, shift) do x1, x2
    x1[1] * x2[2] - x2[1] * x1[2] 
  end
  area = 0.5 * abs(sum(area_components))
  return area
end

function get_facet_diameter(p::Polytope{D}, face::Int) where D
  if D == 3
    @notimplemented
  end
  dim = get_dimranges(p)[face+1]
  X = get_face_coordinates(p)[dim]
  if face == 1
    h = map(X) do x
      norm(x[1]-x[2])
    end
  elseif face == 2
    h = 0.0  
    n_sides = length(X...)
    for i in 1:(n_sides-1)
      for j in (i+1):n_sides
        h = max(h, norm(X[1][i] - X[1][j]))
      end
    end
  end
  return h
end

################## 

function FESpaces.renumber_free_and_dirichlet_dof_ids(
  space::FESpaces.PolytopalFESpace,free_dof_ids,dir_dof_ids
)
  @assert num_free_dofs(space) == length(free_dof_ids)
  @assert num_dirichlet_dofs(space) == length(dir_dof_ids)

  k = PosNegReindex(free_dof_ids,dir_dof_ids)
  cell_dof_ids = lazy_map(Broadcasting(k),get_cell_dof_ids(space))

  ndofs = length(free_dof_ids) + length(dir_dof_ids)
  nfree = sum(i -> i > 0, free_dof_ids; init= 0) + sum(i -> i > 0, dir_dof_ids; init= 0)
  ndirichlet = ndofs - nfree

  cell_is_dirichlet = map(I -> any(i -> i < 0, I), cell_dof_ids)
  dirichlet_cells = collect(Int32,findall(cell_is_dirichlet))

  return PolytopalFESpace(
    space.vector_type,
    nfree,
    ndirichlet,
    cell_dof_ids,
    space.fe_basis,
    cell_is_dirichlet,
    space.dirichlet_dof_tag,
    dirichlet_cells,
    space.ntags,
    space.order,
    space.metadata
  )
end
