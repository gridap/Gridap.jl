"""

The exported names are
$(EXPORTS)
"""
module CellData

using Test
using DocStringExtensions
using FillArrays

using Gridap.Helpers
using Gridap.Arrays
using Gridap.Algebra
using Gridap.TensorValues
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Integration

import Gridap.Arrays: lazy_append
import Gridap.Arrays: get_array
import Gridap.Arrays: evaluate!
import Gridap.Fields: gradient
import Gridap.Fields: ∇∇
import Gridap.Fields: integrate
import Gridap.Fields: grad2curl
import Gridap.Geometry: num_cells
import Gridap.Geometry: get_triangulation

import Gridap.TensorValues: inner, outer, double_contraction, symmetric_part
import LinearAlgebra: det, tr, cross, dot
import Base: inv, abs, abs2, *, +, -, /, adjoint, transpose

export gradient, ∇
export ∇∇
export inner, ⊙, outer, ⊗, double_contraction, ⋅¹, ⋅², symmetric_part
export det, tr, cross, ×, dot, ⋅

export DomainStyle
export ReferenceDomain
export PhysicalDomain
export CellDatum
export get_cell_data
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

export DomainContribution
export num_domains
export get_domains
export get_contribution
export add_contribution!
export LebesgueMeasure

export attach_dirichlet
export attach_constraints_rows
export attach_constraints_cols
export identity_constraints

export CellState
export update!

include("CellDataInterface.jl")

include("CellFields.jl")

include("CellQuadratures.jl")

include("CellStates.jl")

include("DomainContributions.jl")

include("CellDofs.jl")

include("AttachDirichlet.jl")

include("AttachConstraints.jl")

end # module
