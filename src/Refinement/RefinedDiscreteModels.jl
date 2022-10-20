

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
  The hierarchy is stored, allowing for the transfer of dofs 
  between `FESpaces` defined on this model and its parent using a 
  `RefinementTransferoperator`. 
"""
struct RefinedDiscreteModel{Dc,Dp,A<:DiscreteModel{Dc,Dp},B<:DiscreteModel{Dc,Dp},C<:RefinementGlue} <: DiscreteModel{Dc,Dp}
  model  ::A
  parent ::B
  glue   ::C
end

# DiscreteModel API
Geometry.get_grid(model::RefinedDiscreteModel)          = get_grid(model.model)
Geometry.get_grid_topology(model::RefinedDiscreteModel) = get_grid_topology(model.model)
Geometry.get_face_labeling(model::RefinedDiscreteModel) = get_face_labeling(model.model)

# Other getters
get_model(model::RefinedDiscreteModel)  = model.model
get_parent(model::RefinedDiscreteModel) = model.parent
get_glue(model::RefinedDiscreteModel)   = model.glue

function refine(model::DiscreteModel) :: RefinedDiscreteModel
  @abstractmethod
end


# Cartesian builder 

function RefinedCartesianDiscreteModel(domain::Tuple,nC::Int,ref::Int)
  nF = ref*nC
  # Models
  parent = CartesianDiscreteModel(domain,(nC,nC))
  child  = CartesianDiscreteModel(domain,(nF,nF))

  # Glue 
  faces_map      = [Int[],Int[],_create_f2c_cell_map(nC,ref)]
  fcell_child_id = _create_child_map(nC,ref)
  reffe          = LagrangianRefFE(Float64,QUAD,1)
  ref_cell_map   = get_f2c_ref_cell_map(reffe,ref)
  glue = RefinementGlue(faces_map,fcell_child_id,ref_cell_map)

  # RefinedModel
  model = RefinedDiscreteModel(child,parent,glue)
  return model
end

function _create_f2c_cell_map(nC::Int,ref::Int)
  nF = nC*ref

  idx = Tuple.(CartesianIndices((nF,nF)))
  a = map((i,j)->(1+(i-1)÷ref,1+(j-1)÷ref),first.(idx),last.(idx))
  b = map((i,j)->(i-1)*nC+j,first.(a),last.(a))
  return Array(reshape(transpose(b),nF*nF))
end

function _create_child_map(nC::Int,ref::Int)
  nF = nC*ref
  elem = reshape(collect(1:ref*ref),(ref,ref))
  slice = transpose(repeat(elem,nC))
  mat = repeat(slice,nC)
  return Array(reshape(transpose(mat),nF*nF))
end

