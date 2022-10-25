

"""
Refinement glue between two nested triangulations
- `f2c_faces_map`     : Given a fine face gid, returns the gid of the coarse face containing it. 
- `fcell_to_child_id` : Given a fine cell gid, returns the local child id within the coarse cell containing it.
- `f2c_ref_cell_map`  : Ref coordinate map between the fcells and their parent ccell, i.e 
                          for each fine cell defines Φ st. x_c = Φ(x_f)
"""
struct RefinementGlue{A,B,C} <: GridapType
  f2c_faces_map::A
  fcell_to_child_id::B
  f2c_ref_cell_map::C

  function RefinementGlue(f2c_faces_map,
                          fcell_to_child_id,
                          f2c_ref_cell_map)
    A = typeof(f2c_faces_map)
    B = typeof(fcell_to_child_id)
    C = typeof(f2c_ref_cell_map)
    new{A,B,C}(f2c_faces_map,fcell_to_child_id,f2c_ref_cell_map)
  end
end

"""
Ref coordinate map from fcell ref coords to ccell ref coords.
Size -> Number of children per coarse cell.
"""
function get_f2c_ref_cell_map(reffe,ref)
  ref_grid = Geometry.UnstructuredGrid(Visualization.compute_reference_grid(reffe,ref))
  return get_cell_map(ref_grid)
end

"""
Ref coordinate map between fine and coarse cells. 
Size -> Number of fine cells. 
"""
function get_f2c_ref_coordinate_map(g::RefinementGlue)
  m = Reindex(g.f2c_ref_cell_map)
  return lazy_map(m,g.fcell_to_child_id)
end

#################################################################

"""
  RefinedDiscreteModel

  `DiscreteModel` created by refining another `DiscreteModel`. 
  The refinement hierarchy can be traced backwards by following the 
  `parent` pointer chain. This allows the transfer of dofs 
  between `FESpaces` defined on this model and its ancestors using a 
  `RefinementTransferOperator`.

"""
struct RefinedDiscreteModel{Dc,Dp,A<:DiscreteModel{Dc,Dp},B<:DiscreteModel{Dc,Dp},C<:RefinementGlue} <: DiscreteModel{Dc,Dp}
  model  ::A
  parent ::B
  glue   ::C

  function RefinedDiscreteModel(model::DiscreteModel{Dc,Dp},parent,glue) where {Dc,Dp}
    @check !isa(model,RefinedDiscreteModel)
    A = typeof(model)
    B = typeof(parent)
    C = typeof(glue)
    return new{Dc,Dp,A,B,C}(model,parent,glue)
  end
end

# DiscreteModel API
Geometry.get_grid(model::RefinedDiscreteModel)          = get_grid(model.model)
Geometry.get_grid_topology(model::RefinedDiscreteModel) = get_grid_topology(model.model)
Geometry.get_face_labeling(model::RefinedDiscreteModel) = get_face_labeling(model.model)

# Other getters
get_model(model::RefinedDiscreteModel)  = model.model
get_parent(model::RefinedDiscreteModel) = model.parent
get_parent(model::RefinedDiscreteModel{A,<:RefinedDiscreteModel,C}) where {A,C} = get_model(model.parent)
get_refinement_glue(model::RefinedDiscreteModel) = model.glue

function refine(model::DiscreteModel,args...;kwargs...) :: RefinedDiscreteModel
  @abstractmethod
end

function refine(model::RefinedDiscreteModel,args...;kwargs...)
  ref_model = refine(model.model,args...;kwargs...)
  return RefinedDiscreteModel(ref_model.model,model,ref_model.glue)
end

# Cartesian Mesh refining

function refine(model::CartesianDiscreteModel; num_refinements::Int=2)
  @check num_refinements >= 2
  ref  = num_refinements
  desc = Geometry.get_cartesian_descriptor(model)
  nC   = desc.partition

  @notimplementedif length(nC) != 2
  @notimplementedif any(map(nCi -> nCi != nC[1],nC))

  domain    = _get_cartesian_domain(desc)
  model_ref = CartesianDiscreteModel(domain,ref.*nC)

  # Glue
  faces_map      = [Int[],Int[],_create_f2c_cell_map(nC,ref)]
  fcell_child_id = _create_child_map(nC,ref)
  reffe          = LagrangianRefFE(Float64,QUAD,1)
  ref_cell_map   = get_f2c_ref_cell_map(reffe,ref)
  glue = RefinementGlue(faces_map,fcell_child_id,ref_cell_map)

  # RefinedModel
  return RefinedDiscreteModel(model_ref,model,glue)
end

function _get_cartesian_domain(desc::CartesianDescriptor{D}) where D
  origin = desc.origin
  corner = origin + VectorValue(desc.sizes .* desc.partition)
  domain = Vector{Int}(undef,2*D)
  for d in 1:D
    domain[d*2-1] = origin[d]
    domain[d*2]   = corner[d]
  end
  return Tuple(domain)
end

function _create_f2c_cell_map(nC::Tuple,ref::Int)
  D  = length(nC)
  nF = nC[1]*ref

  idx = Tuple.(CartesianIndices((nF,nF)))
  a = map((i,j)->(1+(i-1)÷ref,1+(j-1)÷ref),first.(idx),last.(idx))
  b = map((i,j)->(i-1)*nC[1]+j,first.(a),last.(a))
  return Array(reshape(transpose(b),nF*nF))
end

function _create_child_map(nC::Tuple,ref::Int)
  nF = nC[1]*ref
  elem = reshape(collect(1:ref*ref),(ref,ref))
  slice = transpose(repeat(elem,nC[1]))
  mat = repeat(slice,nC[1])
  return Array(reshape(transpose(mat),nF*nF))
end

