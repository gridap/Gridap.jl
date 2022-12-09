
struct LinearizedFESpace{M,RM,LF,G} <: SingleFieldFESpace
  model::M
  refined_model::RM
  linear_fe_space::LF
  ref_dofs_grid::G
end

function _linearize_reffe(reffe::Tuple{<:Lagrangian,Any,Any})
   (reffe[1],(reffe[2][1],1),reffe[3])
end 

function _num_cells_x_dim(reffe::Tuple{<:Lagrangian,Any,Any})
   reffe[2][2]
end 

function LinearizedFESpace(model::CartesianDiscreteModel{Dc},
                           reffe::Tuple{<:Lagrangian,Any,Any}; kwargs...) where Dc 
  reffe_linearized=_linearize_reffe(reffe)
  ncells_x_dim=_num_cells_x_dim(reffe)
  println("xxx: $(ncells_x_dim)")
  ref_model=Adaptivity.refine(model,ncells_x_dim)
  linear_fe_space=FESpace(ref_model,reffe_linearized;kwargs...)
  polytopes=get_polytopes(model)
  @assert length(polytopes) == 1
  partition=Tuple([ncells_x_dim for i=1:Dc])
  ref_dofs_grid = UnstructuredGrid(compute_reference_grid(polytopes[1],partition))
  LinearizedFESpace(model,ref_model,linear_fe_space,ref_dofs_grid)
end

# FESpace interface
FESpaces.get_free_dof_ids(f::LinearizedFESpace) = get_free_dof_ids(f.linear_fe_space)
FESpaces.zero_free_values(f::LinearizedFESpace) = zero_free_values(f.linear_fe_space)

struct GlueChildrenMap{T,A} <: Fields.Map
  ref_dofs_grid::A
end 

function Fields.return_cache(m::GlueChildrenMap{T},grouped_dofs_vals::AbstractArray) where T
  refs_dofs_cell_dof_ids=get_cell_node_ids(m.ref_dofs_grid)
  nnodes=num_nodes(m.ref_dofs_grid)
  array=Vector{T}(undef,nnodes)
  CachedArray(array), refs_dofs_cell_dof_ids
end 

function Fields.evaluate!(cache,m::GlueChildrenMap{T},grouped_dofs_vals::AbstractArray) where T
  ca,refs_dofs_cell_dof_ids=cache
  result=ca.array
  @assert length(refs_dofs_cell_dof_ids)==length(grouped_dofs_vals)
  for i=1:length(refs_dofs_cell_dof_ids)
    dof_ids_current_cell=refs_dofs_cell_dof_ids[i]
    vals_current_cell=grouped_dofs_vals[i]
    for (j,e) in enumerate(dof_ids_current_cell)
        result[e]=vals_current_cell[j]
    end 
  end 
  result
end 

function Fields.return_cache(m::GlueChildrenMap,grouped_dofs_vals::AbstractArray{<:ReferenceFEs.LagrangianDofBasis})
  Geometry.get_node_coordinates(m.ref_dofs_grid)
end 

function Fields.evaluate!(cache,m::GlueChildrenMap,grouped_dofs_vals::AbstractArray{<:ReferenceFEs.LagrangianDofBasis})
  nodes_coordinates=cache
  nnodes=length(nodes_coordinates)
  dof_to_node = [ i for i=1:nnodes ]
  dof_to_comp = [ 1 for i=1:nnodes ]
  node_and_comp_to_dof=[ [[i]] for i=1:nnodes ]
  LagrangianDofBasis(nodes_coordinates,
                     dof_to_node,
                     dof_to_comp,
                     node_and_comp_to_dof)
end 

function FESpaces.get_cell_dof_ids(f::LinearizedFESpace)
  glue = f.refined_model.glue
  fine_cell_dof_ids = get_cell_dof_ids(f.linear_fe_space)
  T=Int32
  m=GlueChildrenMap{T,typeof(f.ref_dofs_grid)}(f.ref_dofs_grid)
  fine_cell_dof_ids_glued = lazy_map(m,Adaptivity.f2c_reindex(fine_cell_dof_ids,glue))
end

function FESpaces.get_fe_dof_basis(f::LinearizedFESpace) 
  glue             = f.refined_model.glue
  node_coordinates = Geometry.get_node_coordinates(f.ref_dofs_grid)
  nnodes           = length(node_coordinates)
  dof_to_node = [ i for i=1:nnodes ]
  dof_to_comp = [ 1 for i=1:nnodes ]
  node_and_comp_to_dof=[ i for i=1:nnodes ]
  dof_basis=LagrangianDofBasis(node_coordinates,
                               dof_to_node,
                               dof_to_comp,
                               node_and_comp_to_dof)
  trian=get_triangulation(f)
  CellDof(Fill(dof_basis,num_cells(trian)),
          trian,
          ReferenceDomain())
end 

Geometry.get_triangulation(f::LinearizedFESpace) = Triangulation(f.model)
FESpaces.get_dof_value_type(f::LinearizedFESpace) = get_dof_value_type(f.linear_fe_space)
FESpaces.get_vector_type(f::LinearizedFESpace) = get_vector_type(f.linear_fe_space)
FESpaces.get_cell_is_dirichlet(f::LinearizedFESpace) = @notimplemented
FESpaces.ConstraintStyle(::Type{<:LinearizedFESpace}) = FESpaces.UnConstrained()

# SingleFieldFESpace interface
FESpaces.get_dirichlet_dof_ids(f::LinearizedFESpace) = get_dirichlet_dof_ids(f.linear_fe_space)
FESpaces.num_dirichlet_tags(f::LinearizedFESpace) = num_dirichlet_tags(f.linear_fe_space)
FESpaces.zero_dirichlet_values(f::LinearizedFESpace) = zero_dirichlet_values(f.linear_fe_space)
FESpaces.get_dirichlet_dof_tag(f::LinearizedFESpace) = get_dirichlet_dof_tag(f.linear_fe_space)

function FESpaces.scatter_free_and_dirichlet_values(f::LinearizedFESpace,free_values,dirichlet_values)
  @check eltype(free_values) == eltype(dirichlet_values) """\n
  The entries stored in free_values and dirichlet_values should be of the same type.

  This error shows up e.g. when trying to build a FEFunction from a vector of integers
  if the Dirichlet values of the underlying space are of type Float64, or when the
  given free values are Float64 and the Dirichlet values ComplexF64.
  """
  cell_dof_ids = get_cell_dof_ids(f)
  lazy_map(Broadcasting(PosNegReindex(free_values,dirichlet_values)),cell_dof_ids)
end

function FESpaces.gather_free_and_dirichlet_values!(free_vals,dirichlet_vals,f::LinearizedFESpace,cell_vals)
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

function gather_dirichlet_values!(dirichlet_vals,f::LinearizedFESpace,cell_vals)
  cell_dofs = get_cell_dof_ids(f)
  cache_vals = array_cache(cell_vals)
  cache_dofs = array_cache(cell_dofs)
  free_vals = zero_free_values(f)
  cells = f.dirichlet_cells
  _free_and_dirichlet_values_fill!(
    free_vals,
    dirichlet_vals,
    cache_vals,
    cache_dofs,
    cell_vals,
    cell_dofs,
    cells)
  dirichlet_vals
end

function  _free_and_dirichlet_values_fill!(
  free_vals,
  dirichlet_vals,
  cache_vals,
  cache_dofs,
  cell_vals,
  cell_dofs,
  cells)

  for cell in cells
    vals = getindex!(cache_vals,cell_vals,cell)
    dofs = getindex!(cache_dofs,cell_dofs,cell)
    for (i,dof) in enumerate(dofs)
      val = vals[i]
      if dof > 0
        free_vals[dof] = val
      elseif dof < 0
        dirichlet_vals[-dof] = val
      else
        @unreachable "dof ids either positive or negative, not zero"
      end
    end
  end
end

function FESpaces.get_fe_basis(f::LinearizedFESpace)
  fe_basis = get_fe_basis(f.linear_fe_space)
  num_cells=length(f.refined_model.glue.refinement_rules)
  cell_array=lazy_map(FineToCoarseBasis,
                      Adaptivity.f2c_reindex(fe_basis,f.refined_model.glue),
                      Fill(f.ref_dofs_grid,num_cells),
                      f.refined_model.glue.refinement_rules)
  FESpaces.SingleFieldFEBasis(cell_array,get_triangulation(f),FESpaces.TestBasis(),ReferenceDomain())
end

struct FineToCoarseBasis{A<:AbstractArray{<:AbstractArray{<:Field}},
                         B,
                         C<:Adaptivity.RefinementRule} <: AbstractVector{Field}
  fine_fields    :: A
  ref_dofs_grid  :: B
  rrule          :: C 
  function FineToCoarseBasis(fine_fields,ref_dofs_grid,rrule::Adaptivity.RefinementRule)
    @check length(fine_fields)   == Adaptivity.num_subcells(rrule)
    @check num_cells(ref_dofs_grid) == Adaptivity.num_subcells(rrule)
    A = typeof(fine_fields)
    B = typeof(ref_dofs_grid)
    C = typeof(rrule)
    new{A,B,C}(fine_fields,ref_dofs_grid,rrule)
  end
end

function Base.length(f::FineToCoarseBasis)
  Geometry.num_nodes(f.ref_dofs_grid)
end 
function Base.size(f::FineToCoarseBasis)
  (Geometry.num_nodes(f.ref_dofs_grid),)
end 
function Base.getindex(f::FineToCoarseBasis,i)
  ConstantField(0.0)  
end

function Geometry.return_cache(a::FineToCoarseBasis,x::AbstractArray{<:Point})
  fields, x_to_cell = a.fine_fields, a.rrule.x_to_cell
  cmaps = Adaptivity.get_inverse_cell_map(a.rrule)

  xi_cache = array_cache(x)
  fi_cache = array_cache(fields)
  mi_cache = array_cache(cmaps)

  xi = getindex!(xi_cache,x,1)
  child_id = x_to_cell(xi)
  mi = getindex!(mi_cache,cmaps,child_id)
  fi = getindex!(fi_cache,fields,child_id)

  zi_cache = Fields.return_cache(mi,xi)
  zi = evaluate!(zi_cache,mi,xi)

  yi_type  = Fields.return_type(fi,zi)
  yi_cache = Fields.return_cache(fi,zi)
  y_cache  = Arrays.CachedArray(eltype(yi_type),2)

  return fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache
end

function Geometry.evaluate!(cache,a::FineToCoarseBasis{<:AbstractArray{<:AbstractArray{<:Field}}},x::AbstractArray{<:Point})
  fi_cache, mi_cache, xi_cache, zi_cache, yi_cache, y_cache = cache
  fields, x_to_cell = a.fine_fields, a.rrule.x_to_cell
  cmaps = get_inverse_cell_map(a.rrule)

  Arrays.setsize!(y_cache, (length(x),Geometry.num_nodes(a.ref_dofs_grid)))
  
  y_cache.array .= 0.0 

  cell_node_ids = Geometry.get_cell_node_ids(a.ref_dofs_grid)
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
