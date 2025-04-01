
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

function PolytopalFESpace(
  t::Triangulation,args...;trian=nothing,kwargs...
)
  @check isnothing(trian)
  model = get_active_model(t)
  PolytopalFESpace(model,args...;trian=t,kwargs...)
end

function PolytopalFESpace(
  model::DiscreteModel,::Type{T},order::Integer;
  space=:P,vector_type=nothing,
  kwargs...
) where {T}
  D = num_cell_dims(model)
  prebasis = Polynomials.MonomialBasis{D}(T, order, _filter_from_space(space))
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
  dirichlet_tags=Int[],
  dirichlet_masks=nothing
)
  Dc = num_cell_dims(model)

  cell_polytopes = Geometry.get_cell_polytopes(model)
  if eltype(cell_polytopes) <: ReferenceFEs.GeneralPolytope
    cell_change = lazy_map(centroid_map, cell_polytopes)
    cell_shapefuns = lazy_map(Broadcasting(∘), cell_prebasis, cell_change)
    domain_style = PhysicalDomain()
  else
    cell_shapefuns = cell_prebasis
    domain_style = ReferenceDomain()
  end
  fe_basis = SingleFieldFEBasis(cell_shapefuns, trian, TestBasis(), domain_style)
  
  order = maximum(Polynomials.get_order,cell_prebasis)
  ncomps = num_components(return_type(first(cell_prebasis)))
  if !isnothing(dirichlet_masks)
    @check length(dirichlet_masks) == ncomps
    dirichlet_components = dirichlet_masks
  else
    dirichlet_components = fill(false,ncomps)
  end

  ntags = length(dirichlet_tags)
  cell_to_tag = get_face_tag_index(labels,dirichlet_tags,Dc)
  cell_is_dirichlet = map(!isequal(UNSET),cell_to_tag)
  ctype_to_prebasis, cell_to_ctype = compress_cell_data(cell_prebasis)
  ctype_to_conformity = map(MonomialDofConformity,ctype_to_prebasis)
  ctype_to_ldof_to_comp = map(c -> c.dof_to_comp,ctype_to_conformity)
  cell_dof_ids, nfree, ndir, dirichlet_dof_tag, dirichlet_cells = compute_discontinuous_cell_dofs(
    cell_to_ctype, ctype_to_ldof_to_comp, cell_to_tag, dirichlet_components
  )

  metadata = ctype_to_conformity
  return PolytopalFESpace(
    vector_type,nfree,ndir,cell_dof_ids,fe_basis,
    cell_is_dirichlet,dirichlet_dof_tag,dirichlet_cells,ntags,
    order,metadata
  )
end

struct MonomialDofConformity{D,V}
  orders::NTuple{D,Int}
  terms::Vector{CartesianIndex{D}}
  dof_to_term::Vector{Int}
  dof_to_comp::Vector{Int}
  term_and_comp_to_dof::Vector{V}
end

function MonomialDofConformity(basis::Polynomials.MonomialBasis)
  T = return_type(basis)
  nterms = length(basis.terms)
  dof_to_term, dof_to_comp, term_and_comp_to_dof = ReferenceFEs._generate_dof_layout_node_major(T,nterms)
  MonomialDofConformity(
    basis.orders, basis.terms, dof_to_term, dof_to_comp, term_and_comp_to_dof
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

function get_facet_measure(p::Polytope{D}, face::Int) where D
  measures = Float64[]
  if D == 3
    @notimplemented
  elseif isa(p, ExtrusionPolytope{2})
    if p == QUAD 
      perm = [1,2,4,3]
    elseif p == TRI
      perm = [1,2,3]
    end
  elseif isa(p, Polygon)   
    perm = collect(1:length(p.edge_vertex_graph))
  end

  dim = get_dimranges(p)[face+1]
  face_ents = get_face_coordinates(p)[dim]
  if face == 0
    for entity in face_ents
      push!(measures, 0.0)
    end
  elseif face == 1
    for entity in face_ents
      p1, p2 = entity
      push!(measures, norm(p2-p1))
    end
  elseif face == 2
    face_ents = map(Reindex(face_ents...),perm)
    area = shoelace(face_ents)
    push!(measures, area)
  end
  return measures
end

function get_facet_centroid(p::Polytope{D}, face::Int) where D

  if D == 3
    @notimplemented
  end

  dim = get_dimranges(p)[face+1]
  face_coords = get_face_coordinates(p)[dim]
  if isa(p, ExtrusionPolytope{2}) || isa(p, ExtrusionPolytope{1})
    if face == 1 || face == 2
      centroid = mean.(face_coords)
    end
  elseif isa(p, Polygon)
    perm = collect(1:length(p.edge_vertex_graph))
    if face == 1
      centroid = mean.(face_coords)
    elseif face == 2
      ents = map(Reindex(face_coords...),perm)
      shift = circshift(ents, -1)

      components_x = map(ents, shift) do x1, x2
        ( x1[1] + x2[1] ) * ( x1[1] * x2[2] - x2[1] * x1[2] )
      end
      components_y = map(ents, shift) do x1, x2
        ( x1[2] + x2[2] ) * ( x1[1] * x2[2] - x2[1] * x1[2] )
      end
      
      area = get_facet_measure(p, face)
      centroid_x = (1 ./ (6*area)) * sum(components_x)
      centroid_y = (1 ./ (6*area)) * sum(components_y)        
      centroid = VectorValue{2, Float64}(centroid_x..., centroid_y...)
    end
  end
  return centroid
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
