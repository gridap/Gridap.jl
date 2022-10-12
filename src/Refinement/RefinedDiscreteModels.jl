

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
function get_f2c_ref_cell_map(reffe)
  ref_grid = Geometry.UnstructuredGrid(Visualization.compute_reference_grid(reffe,2))
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

struct RefinedDiscreteModel{Dc,Dp,A,B} <: DiscreteModel{Dc,Dp}
  model  ::A
  parent ::B
  glue   ::RefinementGlue

  function RefinedDiscreteModel(model::DiscreteModel{Dc,Dp},parent::DiscreteModel{Dc,Dp},glue) where {Dc,Dp}
    A = typeof(model)
    B = typeof(parent)
    return new{Dc,Dp,A,B}(model,parent,glue)
  end
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
