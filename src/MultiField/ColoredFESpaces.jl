
"""
    abstract type ColoredFESpace <: SingleFieldFESpace end

  FESpace with a color for each dof. At assembly time, each color gets assembled 
  into a separate block of the global matrix.
"""
abstract type ColoredFESpace <: SingleFieldFESpace end

"""
    num_colors(W::ColoredFESpace)

  Returns the number of colors in the `ColoredFESpace`.
"""
function num_colors(::ColoredFESpace)
  @abstractmethod
end

"""
    num_dofs_per_color(W::ColoredFESpace)

  Returns the number of dofs for each color in the `ColoredFESpace`, as an array.
"""
function num_dofs_per_color(::ColoredFESpace)
  @abstractmethod
end

"""
    get_dof_color(W::ColoredFESpace)

  Returns the color of each dof in the `ColoredFESpace`.
"""
function get_dof_color(::ColoredFESpace)
  @abstractmethod
end

"""
    get_cell_dof_color(W::ColoredFESpace)

  Returns the color of each dof in the `ColoredFESpace`, as a cell-wise `Table`.
"""
function get_cell_dof_color(::ColoredFESpace)
  @abstractmethod
end

function get_cell_dof_color(W::ColoredFESpace,trian::Triangulation)
  FESpaces.get_cell_fe_data(get_cell_dof_color,W,trian)
end

"""
    get_dof_to_colored_dof(W::ColoredFESpace)

  Returns the id of each dof within its color group.
"""
function get_dof_to_colored_dof(W::ColoredFESpace)
  get_dof_to_colored_dof(num_dofs_per_color(W),get_dof_color(W))
end

function get_dof_to_colored_dof(color_to_ndofs,dof_to_color)
  dof_to_colored_dof = Vector{Int32}(undef,length(dof_to_color))
  color_counters = fill(0,length(color_to_ndofs))
  for (dof,color) in enumerate(dof_to_color)
    if color > 0
      color_counters[color] += 1
      dof_to_colored_dof[dof] = color_counters[color]
    else
      dof_to_colored_dof[dof] = -1
    end
  end
  @assert color_counters == color_to_ndofs
  return dof_to_colored_dof
end

"""
    struct SingleColorFESpace <: ColoredFESpace end

  A colored FEspace with a single color for all dofs. Mostly used in combination with 
  `MultiColorFESpace`s within `MultiFieldFESpace`s. 
"""
struct SingleColorFESpace{V,A} <: ColoredFESpace
  space :: A
  function SingleColorFESpace(space::SingleFieldFESpace)
    V = get_vector_type(space)
    A = typeof(space)
    new{V,A}(space)
  end
end

# SingleFieldFESpace API

FESpaces.get_fe_basis(W::SingleColorFESpace)       = get_fe_basis(W.space)
FESpaces.get_trial_fe_basis(W::SingleColorFESpace) = get_trial_fe_basis(W.space)
FESpaces.ConstraintStyle(W::SingleColorFESpace)    = ConstraintStyle(W.space)
CellData.get_triangulation(W::SingleColorFESpace)  = get_triangulation(W.space)
FESpaces.get_fe_dof_basis(W::SingleColorFESpace)   = get_fe_dof_basis(W.space)
FESpaces.get_free_dof_ids(W::SingleColorFESpace)   = get_free_dof_ids(W.space)
FESpaces.ConstraintStyle(::Type{<:SingleColorFESpace{V,A}}) where {V,A} = ConstraintStyle(A)
FESpaces.get_vector_type(::SingleColorFESpace{V}) where V = V

function FESpaces.get_cell_dof_ids(W::SingleColorFESpace)
  cell_dof_ids    = get_cell_dof_ids(W.space)
  cell_dof_colors = get_cell_dof_colors(W)
  return lazy_map(ColorMap(),cell_dof_ids,cell_dof_colors)
end

function FESpaces.get_cell_dof_ids(W::SingleColorFESpace,ttrian::Triangulation)
  cell_dof_ids    = get_cell_dof_ids(W.space,ttrian)
  cell_dof_colors = get_cell_dof_color(W,ttrian)
  return lazy_map(ColorMap(),cell_dof_ids,cell_dof_colors)
end

# ColoredFESpace API

num_colors(W::SingleColorFESpace) = 1
num_dofs_per_color(W::SingleColorFESpace) = [num_free_dofs(W.space)]

function get_dof_color(W::SingleColorFESpace)
  return Fill(one(Int8),num_free_dofs(W.space))
end

function get_cell_dof_color(W::SingleColorFESpace)
  cell_dof_ids = get_cell_dof_ids(W.space)
  return Table(get_dof_color(W),cell_dof_ids.ptrs)
end

# This is bad, but it is the only way to make it work for now with MultiConstantFESpaces
function get_cell_dof_color(W::SingleColorFESpace,ttrian::Triangulation)
  cell_dof_ids = get_cell_dof_ids(W.space,ttrian)
  return Table(get_dof_color(W),cell_dof_ids.ptrs)
end

function get_dof_to_colored_dof(W::SingleColorFESpace)
  return Base.OneTo(Int32(num_free_dofs(W.space)))
end


"""
    struct MultiColorFESpace <: ColoredFESpace

    MultiColorFESpace(model,reffe,tags[,labels]; kwargs...)

  A `FESpace` whose dofs are split into disjoint sets, each one associated to a domain tag.
  Free dofs which do not belong to any tag are rejected, and set as homogeneous dirichlet.

  Contains the following fields:
    - `space`: the underlying `FESpace`
    - `tags`: the tags associated to each set of dofs
    - `tag_to_ndofs`: the number of dofs associated to each tag
    - `dof_to_tag`: the tag associated to each dof
    - `dof_to_pdof`: the position of each dof within its tag group
"""
struct MultiColorFESpace{V,A} <: ColoredFESpace
  space        :: A
  tags         :: Vector{String}
  tag_to_ndofs :: Vector{Int64}
  dof_to_tag   :: Vector{Int8}
  dof_to_pdof  :: Vector{Int32}
  function MultiColorFESpace(
    space::SingleFieldFESpace,
    tags::Vector{String},
    tag_to_ndofs::Vector{Int64},
    dof_to_tag::Vector{Int8},
    dof_to_pdof::Vector{Int32}
  )
    V = typeof(mortar(map(n->zeros(n),tag_to_ndofs)))
    A = typeof(space)
    new{V,A}(space,tags,tag_to_ndofs,dof_to_tag,dof_to_pdof)
  end
end

function MultiColorFESpace(
  model ::DiscreteModel,
  reffe ::Tuple{<:Gridap.FESpaces.ReferenceFEName,Any,Any},
  tags  ::Vector{String},
  labels::FaceLabeling = get_face_labeling(model);
  kwargs...
)
  space = FESpace(model,reffe;kwargs...)
  return MultiColorFESpace(space,reffe,tags,labels;kwargs...)
end

function MultiColorFESpace(
  space ::SingleFieldFESpace,
  reffe ::Tuple{<:Gridap.FESpaces.ReferenceFEName,Any,Any},
  tags  ::Vector{String},
  labels::FaceLabeling;
  kwargs...
)
  model = get_background_model(get_triangulation(space))
  cell_conf = _cell_conformity(model,reffe;kwargs...)
  tag_to_ndofs, dof_to_tag = get_dof_to_tag(space,cell_conf,tags,labels)
  dof_to_pdof = get_dof_to_colored_dof(tag_to_ndofs,dof_to_tag)
  return MultiColorFESpace(space,tags,tag_to_ndofs,dof_to_tag,dof_to_pdof)
end

# FESpace interface

FESpaces.get_fe_basis(W::MultiColorFESpace)       = get_fe_basis(W.space)
FESpaces.get_trial_fe_basis(W::MultiColorFESpace) = get_trial_fe_basis(W.space)
FESpaces.ConstraintStyle(W::MultiColorFESpace)    = ConstraintStyle(W.space)
CellData.get_triangulation(W::MultiColorFESpace)  = get_triangulation(W.space)
FESpaces.get_fe_dof_basis(W::MultiColorFESpace)   = get_fe_dof_basis(W.space)
FESpaces.ConstraintStyle(::Type{<:MultiColorFESpace{V,A}}) where {V,A} = ConstraintStyle(A)
FESpaces.get_vector_type(::MultiColorFESpace{V}) where V = V

function FESpaces.get_free_dof_ids(W::MultiColorFESpace)
  return blockedrange(W.tag_to_ndofs)
end

function FESpaces.get_cell_dof_ids(W::MultiColorFESpace)
  cell_dof_ids  = get_cell_dof_ids(W.space)
  ndir = num_dirichlet_dofs(W.space)
  if ndir > 0
    data = lazy_map(PosNegReindex(W.dof_to_pdof,collect(Int32,-2:-1:-(ndir+1))),cell_dof_ids.data)
  else
    data = lazy_map(Reindex(W.dof_to_pdof),cell_dof_ids.data)
  end
  cell_pdof_ids = Table(data,cell_dof_ids.ptrs)
  cell_dof_to_tag = get_cell_dof_color(W)
  return lazy_map(ColorMap(),cell_pdof_ids,cell_dof_to_tag)
end

# ColoredFESpace interface

num_colors(W::MultiColorFESpace) = length(W.tags)
num_dofs_per_color(W::MultiColorFESpace) = W.tag_to_ndofs
get_dof_color(W::MultiColorFESpace) = W.dof_to_tag
get_dof_to_colored_dof(W::MultiColorFESpace) = W.dof_to_pdof

function get_cell_dof_color(W::MultiColorFESpace)
  cell_dof_ids = get_cell_dof_ids(W.space)
  dof_colors = get_dof_color(W)
  ndir = num_dirichlet_dofs(W.space)
  if ndir > 0
    data = lazy_map(PosNegReindex(dof_colors,fill(Int8(-1),ndir)),cell_dof_ids.data)
  else
    data = lazy_map(Reindex(dof_colors),cell_dof_ids.data)
  end
  return Table(data,cell_dof_ids.ptrs)
end

function get_color_cell_dof_ids(W::MultiColorFESpace,color::Integer)
  cell_dof_ids = get_cell_dof_ids(W.space)
  dof_colors = get_dof_color(W)
  dof_to_pdof = lazy_map(ColorMask(color),W.dof_to_pdof,dof_colors)

  ndir = num_dirichlet_dofs(W.space)
  if ndir > 0
    data = lazy_map(PosNegReindex(dof_to_pdof,collect(Int32,-2:-1:-(ndir+1))),cell_dof_ids.data)
  else
    data = lazy_map(Reindex(dof_to_pdof),cell_dof_ids.data)
  end
  return Table(data,cell_dof_ids.ptrs)
end

# Assembly

function FESpaces.SparseMatrixAssembler(
  mat,vec,
  trial::MultiColorFESpace,
  test ::MultiColorFESpace,
  strategy::AssemblyStrategy=DefaultAssemblyStrategy()
)
  matrix_builder = SparseMatrixBuilder(mat)
  vector_builder = ArrayBuilder(vec)

  # Count block rows/cols
  block_rows = blocks(get_free_dof_ids(trial))
  block_cols = blocks(get_free_dof_ids(test))
  @assert length(block_rows) == length(block_cols)
  NB = length(block_rows); SB = Tuple(fill(1,NB)); P  = Tuple(collect(1:NB))

  # Create block assemblers
  block_idx = CartesianIndices((NB,NB))
  block_assemblers = map(block_idx) do idx
    rows = block_rows[idx[1]]
    cols = block_cols[idx[2]]
    FESpaces.GenericSparseMatrixAssembler(matrix_builder,vector_builder,rows,cols,strategy)
  end

  return MultiField.BlockSparseMatrixAssembler{NB,NB,SB,P}(block_assemblers)
end

# This should be moved to Gridap/BlockAssemblers
FESpaces.map_cell_rows(strategy::MatrixBlock{FESpaces.DefaultAssemblyStrategy},cell_ids) = cell_ids
FESpaces.map_cell_cols(strategy::MatrixBlock{FESpaces.DefaultAssemblyStrategy},cell_ids) = cell_ids

# Auxiliary functions

function _cell_conformity(
  model::DiscreteModel,
  reffe::Tuple{<:Gridap.FESpaces.ReferenceFEName,Any,Any};
  conformity=nothing, kwargs...
)
  basis, reffe_args, reffe_kwargs = reffe
  cell_reffe = ReferenceFE(model,basis,reffe_args...;reffe_kwargs...)
  conformity = Conformity(Gridap.Arrays.testitem(cell_reffe),conformity)
  return CellConformity(cell_reffe,conformity)
end

function Geometry.get_face_labeling(space::FESpace)
  return get_face_labeling(get_background_model(get_triangulation(space)))
end

function get_dof_to_tag(
  space::FESpace,
  cell_conf::CellConformity,
  tags::Vector{String},
  labels::FaceLabeling = get_face_labeling(space)
)
  ntags = length(tags)
  ndofs = num_free_dofs(space)

  reject_ndofs = ndofs
  tag_to_ndofs = Vector{Int64}(undef,ntags)
  dof_to_tag   = fill(Int8(-1),ndofs)
  for (i,tag) in enumerate(tags)
    dof_to_mask = get_dof_mask(space,cell_conf,tag,labels)
    dof_to_tag[dof_to_mask] .= i
    i_ndofs = sum(dof_to_mask)
    tag_to_ndofs[i] = i_ndofs
    reject_ndofs   -= i_ndofs
  end
  @assert sum(tag_to_ndofs) == (ndofs-reject_ndofs) "There is overlapping between the tags!"
  return tag_to_ndofs, dof_to_tag
end

function get_dof_mask(
  space::FESpace,
  cell_conf::CellConformity,
  tag::String,
  labels::FaceLabeling = get_face_labeling(space)
)
  model  = get_background_model(get_triangulation(space))
  topo   = get_grid_topology(model)

  Dc = num_cell_dims(model)
  cell_dof_ids = get_cell_dof_ids(space)
  d_to_dface_to_mask = [get_face_mask(labels,tag,d) for d in 0:Dc]
  d_to_cell_to_dface = [Geometry.get_faces(topo,Dc,d) for d in 0:Dc]

  cell_dof_ids_cache = array_cache(cell_dof_ids)
  d_to_cell_to_dface_cache = map(array_cache,d_to_cell_to_dface)
  dof_to_mask = fill(false,num_free_dofs(space))
  for (cell,ctype) in enumerate(cell_conf.cell_ctype)
    dofs = getindex!(cell_dof_ids_cache,cell_dof_ids,cell)
    lface_own_ldofs = cell_conf.ctype_lface_own_ldofs[ctype]
    for d in 0:Dc
      offset = cell_conf.d_ctype_offset[d+1][ctype]
      dfaces = getindex!(d_to_cell_to_dface_cache[d+1],d_to_cell_to_dface[d+1],cell)
      for (lface,dface) in enumerate(dfaces)
        ldofs = lface_own_ldofs[offset+lface]
        for ldof in ldofs
          dof = dofs[ldof]
          if dof > 0 # Avoid dirichlet dofs
            dof_to_mask[dof] |= d_to_dface_to_mask[d+1][dface]
          end
        end
      end
    end
  end

  return dof_to_mask
end
