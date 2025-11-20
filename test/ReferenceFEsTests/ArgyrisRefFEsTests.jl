module ArgyrisRefFEsTest

using Test
using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Helpers
using Gridap.FESpaces
using Gridap.TensorValues
using Gridap.CellData
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Fields

using LinearAlgebra

T = Float64

# Argyris DOF tests TODO

# validate the vertices dofs
# polynomials [ x², y², x²y² ]
#x2_y2 = linear_combination(T[ (j,i)∈((1, 3),(2,7),(3,9)) for i in 1:9, j in 1:3], MonomialBasis(Val(2),T,2))
#shapefuns_monom_coordinates = evaluate(dofs, x2_y2)


# RefFE tests

argyris5 = ArgyrisRefFE(T, TRI, 5)
dofs = get_dof_basis(argyris5)
shapefuns = get_shapefuns(argyris5)
@test norm(evaluate(dofs, shapefuns) - I) <= 1.e-14

@test get_face_own_dofs(argyris5) == Array{Int}[
#  φ  ∂x(φ) ∂y(φ) ∂xx(φ) ∂xy(φ) ∂yy(φ) at v1
  [1, 4,    5,    13,    14,    15],
#  φ  ∂x(φ) ∂y(φ) ∂xx(φ) ∂xy(φ) ∂yy(φ) at v2
  [2, 6,    7,    16,    17,    18],
#  φ  ∂x(φ) ∂y(φ) ∂xx(φ) ∂xy(φ) ∂yy(φ) at v3
  [3, 8,    9,    19,    20,    21],
#  ∫_eᵢ ∂ₙφ deᵢ at edges e₁,e₂,e₂
  [10], [11], [12],
  []
]
@test get_face_dofs(argyris5) == Array{Int}[
  [1, 4, 5, 13, 14, 15], [2, 6, 7, 16, 17, 18], [3, 8, 9, 19, 20, 21],
  [1, 4, 5, 13, 14, 15, 2, 6, 7, 16, 17, 18, 10],
  [1, 4, 5, 13, 14, 15, 3, 8, 9, 19, 20, 21, 11],
  [2, 6, 7, 16, 17, 18, 3, 8, 9, 19, 20, 21, 12],
  [1, 4, 5, 13, 14, 15, 2, 6, 7, 16, 17, 18, 3, 8, 9, 19, 20, 21, 10, 11, 12]
]

argyris7 = ArgyrisRefFE(T, TRI, 7)
dofs = get_dof_basis(argyris7)
shapefuns = get_shapefuns(argyris7)
@test norm(evaluate(dofs, shapefuns) - I) <= 1.e-13


V = VectorValue{2,T}
argyris5_q = ArgyrisRefFE(V, TRI, 5)
dofs = get_dof_basis(argyris5_q)
shapefuns = get_shapefuns(argyris5_q)
@test norm(evaluate(dofs, shapefuns) - I) <= 1.e-13

@test get_face_own_dofs(argyris5_q) == Array{Int}[
#  φ₁ φ₂ ∂xφ₁ ∂xφ₁ ∂yφ₁ ∂yφ₂ ∂xxφ₁ ∂xxφ₂ ∂xyφ₁ ∂xyφ₂ ∂yyφ₁ ∂yyφ₂) at v1
  [1, 2, 7,   8,   9,   10,  25,   26,   27,   28,   29,   30,  ],
#  φ₁ φ₂ ∂xφ₁ ∂xφ₁ ∂yφ₁ ∂yφ₂ ∂xxφ₁ ∂xxφ₂ ∂xyφ₁ ∂xyφ₂ ∂yyφ₁ ∂yyφ₂) at v2
  [3, 4, 11,  12,  13,  14,  31,   32,   33,   34,   35,   36   ],
#  φ₁ φ₂ ∂xφ₁ ∂xφ₁ ∂yφ₁ ∂yφ₂ ∂xxφ₁ ∂xxφ₂ ∂xyφ₁ ∂xyφ₂ ∂yyφ₁ ∂yyφ₂) at v3
  [5, 6, 15,  16,  17,  18,  37,   38,   39,   40,   41,   42   ],
#  ∫_eᵢ ∂ₙφ deᵢ at edges e₁,e₂,e₂
  [19, 20], [21, 22], [23, 24],
  []
]

argyris7_q = ArgyrisRefFE(V, TRI, 7)
dofs = get_dof_basis(argyris7_q)
shapefuns = get_shapefuns(argyris7_q)
@test norm(evaluate(dofs, shapefuns) - I) <= 1.e-13

# Physical mapping of the element TODO

# TODO optimize return_cache(σ,φ,μ,dt::FaceMeasure) (24% dof creation time)
# TODO optimize _compute_reffaces_and_face_types(::Polytope) (17% dof creation time)
# TODO optimize FaceMeasur creation (15% dof creation time)
# Or, add methods to cache and re-uses all the things that don't depend on vertices coords
# - reffaces, face_types
# - measures (have to replace the fmaps)
# - face_nodes, face_own_moms
# -
# - V, op_φ

end
