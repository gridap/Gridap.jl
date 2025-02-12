
struct PolytopalFESpace{V,M} <: SingleFieldFESpace
  vector_type::Type{V}
  ndofs::Int
  cell_dofs_ids::AbstractArray
  fe_basis::CellField
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
  trian::Triangulation{Dc},::Type{T},order::Integer; space = :P, vector_type = Vector{T}
) where {Dc,T}
  ncells = num_cells(trian)
  prebasis = Polynomials.MonomialBasis{Dc}(T, order, _filter_from_space(space))
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
  pmin, pmax = get_bounding_box(poly)
  o = VectorValue(ntuple(i->1.0,Val(D)))
  xc = 0.5 * (pmin + pmax)
  # xc = get_facet_centroid(poly, D) # reverting to original
  h = 0.5 * (pmax - pmin)
  return affine_map(TensorValues.diagonal_tensor(o ./ h), -xc ./ h)
end

function Arrays.lazy_map(::typeof(evaluate),a::LazyArray{<:Fill{typeof(centroid_map)}},x::AbstractVector)
  polys = a.args[1]
  lazy_map(CentroidCoordinateChangeMap(),polys,x)
end

function Arrays.evaluate!(cache,::CentroidCoordinateChangeMap,poly::Polytope{D},x::Point{D}) where D
  pmin, pmax = get_bounding_box(poly)
  xc = 0.5 * (pmin + pmax)
  # xc = get_facet_centroid(poly, D) # reverting to original
  h = 0.5 * (pmax - pmin)
  return (x - xc) ./ h
end

function Arrays.return_cache(::CentroidCoordinateChangeMap,poly::Polytope{D},x::AbstractVector{<:Point{D}}) where D
  return CachedArray(similar(x))
end

function Arrays.evaluate!(cache,::CentroidCoordinateChangeMap,poly::Polytope{D},x::AbstractVector{<:Point{D}}) where D
  setsize!(cache,size(x))
  y = cache.array
  pmin, pmax = get_bounding_box(poly)
  xc = 0.5 * (pmin + pmax)
  # xc = get_facet_centroid(poly, D) # reverting to original
  h = 0.5 * (pmax - pmin)
  for i in eachindex(x)
    y[i] = (x[i] - xc) ./ h
  end
  return y
end

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
