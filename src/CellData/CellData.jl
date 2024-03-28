"""

The exported names are
$(EXPORTS)
"""
module CellData

using Test
using DocStringExtensions
using FillArrays

using Gridap.Helpers
using Gridap.Algebra
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry

using NearestNeighbors
using StaticArrays

import Gridap.Arrays: lazy_append
import Gridap.Arrays: get_array
import Gridap.Arrays: evaluate!
import Gridap.Arrays: return_cache
import Gridap.Fields: gradient, DIV
import Gridap.Fields: ∇∇
import Gridap.Fields: integrate
import Gridap.Fields: grad2curl
import Gridap.Geometry: num_cells
import Gridap.Geometry: get_triangulation

import Gridap.TensorValues: inner, outer, double_contraction, symmetric_part
import LinearAlgebra: det, tr, cross, dot, ⋅, rmul!
import Base: inv, abs, abs2, *, +, -, /, adjoint, transpose, real, imag, conj
import Statistics: mean

export gradient, ∇
export ∇∇
export inner, ⊙, outer, ⊗, double_contraction, ⋅¹, ⋅², symmetric_part
export det, tr, cross, ×, dot, ⋅

export DomainStyle
export ReferenceDomain
export PhysicalDomain
export CellDatum
export get_data
export get_triangulation
export change_domain
export test_cell_datum
export CellPoint
export get_cell_points
export CellField
export jump
export mean
export GenericCellField
export CellQuadrature
export Integrand
export ∫
export CellDof
export get_normal_vector
export get_cell_measure
export Interpolable
export KDTreeSearch

export make_inverse_table
export compute_cell_points_from_vector_of_points

export DomainContribution
export num_domains
export get_domains
export get_contribution
export add_contribution!
export Measure

export attach_dirichlet
export attach_constraints_rows
export attach_constraints_cols
export identity_constraints
export get_physical_coordinate

export CellState
export update_state!

export DiracDelta

export SkeletonCellFieldPair 

include("CellDataInterface.jl")

include("CellFields.jl")

include("CellQuadratures.jl")

include("CellStates.jl")

include("DomainContributions.jl")

include("DiracDeltas.jl")

include("CellDofs.jl")

include("AttachDirichlet.jl")

include("AttachConstraints.jl")

include("SkeletonCellFieldPair.jl")

end # module
