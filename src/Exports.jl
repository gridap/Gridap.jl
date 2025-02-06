
macro publish(mod,name)
  quote
    using Gridap.$mod: $name; export $name
  end
end

# Reexport from LinearAlgebra (just for convenience)
using LinearAlgebra:  det, inv, tr, cross, dot, norm, ×, ⋅
export det, inv, tr, cross, dot, norm, ×, ⋅

@publish Helpers GridapType

@publish Algebra solve
@publish Algebra solve!
@publish Algebra symbolic_setup
@publish Algebra numerical_setup
@publish Algebra numerical_setup!
@publish Algebra LUSolver
@publish Algebra BackslashSolver
@publish Algebra zero_initial_guess
@publish Algebra NLSolver
@publish Algebra get_matrix
@publish Algebra get_vector
@publish Algebra jacobian
@publish Algebra hessian

@publish Arrays array_cache
@publish Arrays getindex!
@publish Arrays get_array
@publish Arrays lazy_map
@publish Arrays Reindex
@publish Arrays Broadcasting
@publish Arrays Operation
@publish Arrays print_op_tree
using Gridap.Arrays: ∑; export ∑

@publish TensorValues VectorValue
@publish TensorValues TensorValue
@publish TensorValues inner
@publish TensorValues outer
@publish TensorValues diagonal_tensor
@publish TensorValues num_components
@publish TensorValues num_indep_components
using Gridap.TensorValues: ⊙; export ⊙
using Gridap.TensorValues: ⊗; export ⊗

@publish Fields gradient
@publish Fields ∇
@publish Fields ∇∇
@publish Fields integrate
@publish Fields Point
@publish Fields evaluate
@publish Fields evaluate!
@publish Fields curl
@publish Fields laplacian
@publish Fields divergence
@publish Fields DIV
@publish Fields Δ
@publish Fields ε
@publish Fields symmetric_gradient

@publish ReferenceFEs is_simplex
@publish ReferenceFEs is_n_cube
@publish ReferenceFEs simplexify
@publish ReferenceFEs num_dims
@publish ReferenceFEs num_cell_dims
@publish ReferenceFEs num_point_dims
@publish ReferenceFEs VERTEX
@publish ReferenceFEs SEGMENT
@publish ReferenceFEs TRI
@publish ReferenceFEs QUAD
@publish ReferenceFEs TET
@publish ReferenceFEs HEX
@publish ReferenceFEs WEDGE
@publish ReferenceFEs PYRAMID
@publish ReferenceFEs is_first_order
@publish ReferenceFEs is_Q
@publish ReferenceFEs is_P
@publish ReferenceFEs is_S
@publish ReferenceFEs VERTEX1
@publish ReferenceFEs SEG2
@publish ReferenceFEs TRI3
@publish ReferenceFEs QUAD4
@publish ReferenceFEs TET4
@publish ReferenceFEs HEX8
@publish ReferenceFEs Polytope
@publish ReferenceFEs ReferenceFE
@publish ReferenceFEs Lagrangian
@publish ReferenceFEs RaviartThomas
@publish ReferenceFEs BDM
@publish ReferenceFEs Nedelec
@publish ReferenceFEs ModalC0
@publish ReferenceFEs lagrangian
@publish ReferenceFEs raviart_thomas
@publish ReferenceFEs bdm
@publish ReferenceFEs nedelec
@publish ReferenceFEs modalC0

@publish Geometry get_triangulation
@publish Geometry num_cells
@publish Geometry num_facets
@publish Geometry num_vertices
@publish Geometry num_edges
@publish Geometry num_faces
@publish Geometry Triangulation
@publish Geometry get_cell_coordinates
@publish Geometry get_cell_ref_coordinates
@publish Geometry get_cell_map
@publish Geometry get_glue
@publish Geometry CartesianGrid
@publish Geometry CartesianDiscreteModel
@publish Geometry DiscreteModel
@publish Geometry DiscreteModelFromFile
@publish Geometry num_tags
@publish Geometry num_entities
@publish Geometry get_grid
@publish Geometry get_face_labeling
@publish Geometry add_tag!
@publish Geometry add_tag_from_tags!
@publish Geometry BoundaryTriangulation
@publish Geometry SkeletonTriangulation
@publish Geometry InterfaceTriangulation
@publish Geometry Interior
@publish Geometry Boundary
@publish Geometry Skeleton
@publish Geometry Interface
@publish Geometry move_contributions
@publish Geometry get_background_model
@publish Geometry get_active_model

@publish CellData CellQuadrature
@publish CellData Measure
@publish CellData DomainStyle
@publish CellData ReferenceDomain
@publish CellData PhysicalDomain
@publish CellData get_cell_points
@publish CellData CellField
@publish CellData CellState
@publish CellData jump
@publish CellData mean
@publish CellData update_state!
@publish CellData get_normal_vector
@publish CellData get_tangent_vector
using Gridap.CellData: ∫; export ∫
@publish CellData get_cell_measure
@publish CellData get_physical_coordinate
@publish CellData DiracDelta

@publish FESpaces FESpace
@publish FESpaces TrialFESpace
@publish FESpaces TestFESpace
@publish FESpaces AffineFEOperator
@publish FESpaces LinearFESolver
@publish FESpaces get_free_dof_values
@publish FESpaces get_dirichlet_dof_values
@publish FESpaces num_dirichlet_dofs
@publish FESpaces num_free_dofs
@publish FESpaces num_dirichlet_tags
@publish FESpaces get_free_dof_ids
@publish FESpaces get_dirichlet_dof_ids
@publish FESpaces get_cell_dof_ids
@publish FESpaces get_cell_dof_values
@publish FESpaces get_fe_basis
@publish FESpaces get_trial_fe_basis
@publish FESpaces FEFunction
@publish FESpaces interpolate
@publish FESpaces interpolate_everywhere
@publish FESpaces interpolate_dirichlet
@publish FESpaces assemble_vector
@publish FESpaces assemble_matrix
@publish FESpaces assemble_matrix_and_vector
@publish FESpaces FEOperator
@publish FESpaces FESolver
@publish FESpaces SparseMatrixAssembler
@publish FESpaces FiniteElements
@publish FESpaces ConstantFESpace

@publish MultiField MultiFieldFESpace
@publish MultiField num_fields

@publish Visualization writevtk
@publish Visualization createvtk
@publish Visualization createpvd
@publish Visualization savepvd

@publish ODEs ∂t
@publish ODEs ∂tt
@publish ODEs ForwardEuler
@publish ODEs ThetaMethod
@publish ODEs MidPoint
@publish ODEs BackwardEuler
@publish ODEs GeneralizedAlpha1
@publish ODEs ButcherTableau
@publish ODEs available_tableaus
@publish ODEs RungeKutta
# @publish ODEs GeneralizedAlpha2
# @publish ODEs Newmark
@publish ODEs TransientTrialFESpace
@publish ODEs TransientMultiFieldFESpace
@publish ODEs TransientFEOperator
@publish ODEs TransientIMEXFEOperator
@publish ODEs TransientSemilinearFEOperator
@publish ODEs TransientQuasilinearFEOperator
@publish ODEs TransientLinearFEOperator

# Deprecated / removed

export apply
function apply(args...)
  Helpers.@unreachable """\n
  Function apply has been removed and replaced by lazy_map.
  This error message will be deleted in future versions.
  """
end

export cell_measure
function cell_measure(args...)
  Helpers.@unreachable """\n
  Function cell_measure(a,b) has been removed and replaced by get_cell_measure(a).
  This error message will be deleted in future versions.
  """
end

export restrict
function restrict(args...)
  Helpers.@unreachable """\n
  Function restrict has been removed. The user does not need to explicitly
  restrict to a given Triangulation any more. The code does it undere the hood.
  This error message will be deleted in future versions.
  """
end

export FETerm
function FETerm(args...)
  Helpers.@unreachable """\n
  Function FETerm has been removed. The API for specifying the weak form has changed significantly.
  See the gridap/Tutorials repo for some examples of how to use the new API.
  This error message will be deleted in future versions.
  """
end

export FEEnergy
function FEEnergy(args...)
  Helpers.@unreachable """\n
  Function FEEnergy has been removed. The API for specifying the weak form has changed significantly.
  See the gridap/Tutorials repo for some examples of how to use the new API.
  This error message will be deleted in future versions.
  """
end

export AffineFETerm
function AffineFETerm(args...)
  Helpers.@unreachable """\n
  Function AffineFETerm has been removed. The API for specifying the weak form has changed significantly.
  See the gridap/Tutorials repo for some examples of how to use the new API.
  This error message will be deleted in future versions.
  """
end

export LinearFETerm
function LinearFETerm(args...)
  Helpers.@unreachable """\n
  Function LinearFETerm has been removed. The API for specifying the weak form has changed significantly.
  See the gridap/Tutorials repo for some examples of how to use the new API.
  This error message will be deleted in future versions.
  """
end

export FESource
function FESource(args...)
  Helpers.@unreachable """\n
  Function FESource has been removed. The API for specifying the weak form has changed significantly.
  See the gridap/Tutorials repo for some examples of how to use the new API.
  This error message will be deleted in future versions.
  """
end

@publish FESpaces get_free_values
@publish FESpaces get_dirichlet_values
