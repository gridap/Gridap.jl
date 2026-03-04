
"""
    struct ConstantFESpace <: SingleFieldFESpace
      # private fields
    end

    ConstantFESpace(model::DiscreteModel; vector_type=Vector{Float64}, field_type=Float64)
    ConstantFESpace(trian::Triangulation; vector_type=Vector{Float64}, field_type=Float64)

FESpace that is constant over the provided model/triangulation. Typically used as  
lagrange multipliers. The kwargs `vector_type` and `field_type` are used to specify the
types of the dof-vector and dof-value respectively.
"""
struct ConstantFESpace{V,T,A,B,C,D} <: SingleFieldFESpace
  trian::A
  cell_basis::B
  cell_dof_basis::C
  cell_dof_ids::D

  function ConstantFESpace(
    trian::Triangulation,
    cell_basis::FEBasis,
    cell_dof_basis::CellDof,
    cell_dof_ids::AbstractVector{<:AbstractVector{<:Integer}},
    vector_type::Type{V},
    field_type::Type{T}
  ) where {V,T}
    A, B, C, D = typeof(trian), typeof(cell_basis), typeof(cell_dof_basis), typeof(cell_dof_ids)
    new{V,T,A,B,C,D}(trian, cell_basis, cell_dof_basis, cell_dof_ids)
  end
end

function ConstantFESpace(
  model::DiscreteModel{Dc},
  trian::Triangulation{Dc};
  vector_type::Type{V}=Vector{Float64},
  field_type::Type{T}=Float64
) where {Dc,V,T}
  @assert num_cells(model) == num_cells(trian)

  basis, reffe_args, reffe_kwargs = ReferenceFE(lagrangian,T,0)
  cell_reffe = ReferenceFE(model,basis,reffe_args...;reffe_kwargs...)
  cell_basis_array = lazy_map(get_shapefuns,cell_reffe)

  cell_basis = SingleFieldFEBasis(
    cell_basis_array, trian, TestBasis(), ReferenceDomain()
  )
  cell_dof_basis = CellDof(
    lazy_map(get_dof_basis,cell_reffe),trian,ReferenceDomain()
  )
  cell_dof_ids = Fill(Int32(1):Int32(num_indep_components(field_type)),num_cells(trian))

  return ConstantFESpace(trian, cell_basis, cell_dof_basis, cell_dof_ids, vector_type, field_type)
end

function ConstantFESpace(
  model::Geometry.PolytopalDiscreteModel{Dc},
  trian::Triangulation{Dc};
  vector_type::Type{V}=Vector{Float64},
  field_type::Type{T}=Float64
) where {Dc,V,T}
  @assert num_cells(model) == num_cells(trian)

  ncells = num_cells(trian)
  prebasis = Polynomials.MonomialBasis{Dc}(T, 0)
  cell_basis = SingleFieldFEBasis(
    Fill(prebasis, ncells), trian, TestBasis(), ReferenceDomain()
  )
  cell_dof_basis = CellDof(
    Fill(ReferenceFEs.MockDofBasis([zero(VectorValue{Dc,Float64})]), ncells),trian,ReferenceDomain()
  )
  cell_dof_ids = Fill(Int32(1):Int32(num_indep_components(field_type)),ncells)

  return ConstantFESpace(trian, cell_basis, cell_dof_basis, cell_dof_ids, vector_type, field_type)
end

function ConstantFESpace(model::DiscreteModel; kwargs...)
  trian = Triangulation(model)
  ConstantFESpace(model,trian; kwargs...)
end

function ConstantFESpace(trian::Triangulation; kwargs...)
  model = get_active_model(trian)
  ConstantFESpace(model,trian; kwargs...)
end

TrialFESpace(f::ConstantFESpace) = f

# Delegated functions
get_triangulation(f::ConstantFESpace) = f.trian

ConstraintStyle(::Type{<:ConstantFESpace}) = UnConstrained()

get_dirichlet_dof_values(f::ConstantFESpace{V}) where V = eltype(V)[]

get_fe_basis(f::ConstantFESpace) = f.cell_basis

get_fe_dof_basis(f::ConstantFESpace) = f.cell_dof_basis

get_dof_value_type(f::ConstantFESpace) = Float64

get_free_dof_ids(f::ConstantFESpace) = Base.OneTo(length(f.cell_dof_ids[1]))

get_vector_type(f::ConstantFESpace{V}) where V = V

get_cell_dof_ids(f::ConstantFESpace) = f.cell_dof_ids

get_dirichlet_dof_ids(f::ConstantFESpace) = Base.OneTo(0)

num_dirichlet_tags(f::ConstantFESpace) = 0

get_dirichlet_dof_tag(f::ConstantFESpace) = Int8[]

function scatter_free_and_dirichlet_values(f::ConstantFESpace,fv,dv)
  cell_dof_ids = get_cell_dof_ids(f)
  lazy_map(Broadcasting(PosNegReindex(fv,dv)),cell_dof_ids)
end

function gather_free_and_dirichlet_values!(free_vals, dirichlet_vals, f::ConstantFESpace, cell_vals)
  cell_dofs = get_cell_dof_ids(f)
  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)
  cells = 1:length(cell_vals)

  _free_and_dirichlet_values_fill!(
    free_vals,
    dirichlet_vals,
    cache_vals,
    cache_dofs,
    cell_vals,
    cell_dofs,
    cells)

  (free_vals,dirichlet_vals)
end

############################################################################################

"""
    struct MultiConstantFESpace{N,V,T} <: SingleFieldFESpace end

    MultiConstantFESpace(model::DiscreteModel,tags::Vector,D::Integer)
    MultiConstantFESpace(trians::Vector{<:BoundaryTriangulation{Df}})

Extension of `ConstantFESpace`, representing a FESpace which is constant on each 
of it's N triangulations.
"""
struct MultiConstantFESpace{N,T,V} <: SingleFieldFESpace
  space :: FESpace
  trian :: Triangulation
  face_dofs :: Table{Int32}
end

function MultiConstantFESpace(
  trian::Triangulation, tface_to_subtrian;
  vector_type::Type{V} = Vector{Float64},
  field_type::Type{T} = Float64
) where {T,V}
  N = maximum(tface_to_subtrian)
  reffe = ReferenceFE(lagrangian,T,0)
  space = FESpace(trian,reffe)

  face_dofs = Table(collect(Int32,tface_to_subtrian),Base.OneTo(Int32(num_cells(trian)+1)))
  return MultiConstantFESpace{N,T,V}(space,trian,face_dofs)
end

function MultiConstantFESpace(
  model::DiscreteModel,tags::Vector,D;
  labels = get_face_labeling(model),
  field_type = Float64, vector_type = Vector{Float64}
)
  function onlyonetrue(masks)
    s = sum(masks)
    @assert s <= 1 "Each face must only have one associated tag!"
    return isone(s)
  end

  masks = map(tag -> get_face_mask(labels,tag,D), tags)
  tfaces = findall(map(onlyonetrue,zip(masks...)))
  @check !isempty(tfaces)

  tface_to_subtrian = fill(Int32(0),length(tfaces))
  for (k,mask) in enumerate(masks)
    tface_mask = view(mask,tfaces)
    tface_to_subtrian[tface_mask] .= Int32(k)
  end
  @check !any(iszero,tface_to_subtrian)

  trian = Triangulation(ReferenceFE{D},model,tfaces)
  return MultiConstantFESpace(trian,tface_to_subtrian;field_type,vector_type)
end

function MultiConstantFESpace(
  trians::Vector{<:BoundaryTriangulation{Df}};
  field_type = Float64, vector_type = Vector{Float64}
) where Df
  model = get_background_model(first(trians))
  @check all(get_background_model(trian) === model for trian in trians)

  bgfaces = map(t -> t.glue.face_to_bgface, trians)

  tface_to_mface = bgfaces[1]
  tface_to_subtrian = ones(Int32,length(tface_to_mface))
  for (k,tf) in enumerate(bgfaces[2:end])
    @assert isempty(intersect(tface_to_mface,tf))
    tface_to_mface = vcat(tface_to_mface,tf)
    tface_to_subtrian = vcat(tface_to_subtrian,fill(Int32(k+1),length(tf)))
  end

  trian = Triangulation(ReferenceFE{Df},model,tface_to_mface)
  MultiConstantFESpace(trian,tface_to_subtrian;field_type,vector_type)
end

# SingleFieldFESpace API

get_free_dof_ids(::MultiConstantFESpace{N,T}) where {N,T} = Base.OneTo(Int32(N*num_components(T)))
get_vector_type(::MultiConstantFESpace{N,T,V}) where {N,T,V} = V

ConstraintStyle(::Type{<:MultiConstantFESpace}) = UnConstrained()
get_dirichlet_dof_values(::MultiConstantFESpace{N,T}) where {N,T} = T[]
get_dirichlet_dof_ids(f::MultiConstantFESpace) = Base.OneTo(zero(Int32))
num_dirichlet_tags(f::MultiConstantFESpace) = 0
get_dirichlet_dof_tag(f::MultiConstantFESpace) = Int8[]

get_fe_basis(f::MultiConstantFESpace) = get_fe_basis(f.space)
get_fe_dof_basis(f::MultiConstantFESpace) = get_fe_dof_basis(f.space)

Base.zero(f::MultiConstantFESpace) = zero(f.space)

Geometry.get_triangulation(f::MultiConstantFESpace) = f.trian
get_cell_dof_ids(f::MultiConstantFESpace) = f.face_dofs

function scatter_free_and_dirichlet_values(f::MultiConstantFESpace,fv,dv)
  cell_dof_ids = get_cell_dof_ids(f)
  lazy_map(Broadcasting(PosNegReindex(fv,dv)),cell_dof_ids)
end

function gather_free_and_dirichlet_values!(free_values, dirichlet_values,f::MultiConstantFESpace,cell_vals)
  cell_dofs = get_cell_dof_ids(f)
  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)
  cells = 1:length(cell_vals)

  _free_and_dirichlet_values_fill!(
    free_values,
    dirichlet_values,
    cache_vals,
    cache_dofs,
    cell_vals,
    cell_dofs,
    cells
  )

  (free_values,dirichlet_values)
end
