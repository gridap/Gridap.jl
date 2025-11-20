module HermiteRefFEsTest

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

hermite1 = HermiteRefFE(T, SEGMENT)
dofs = get_dof_basis(hermite1)
shapefuns = get_shapefuns(hermite1)
@test norm(evaluate(dofs, shapefuns) - I) <= 1.e-15

P3_1D_monoms = FEEC_poly_basis(Val(1), T, 3, 0, :P, Monomial)
shapefuns_monom_coordinates = inv(evaluate(dofs, P3_1D_monoms))
@test shapefuns_monom_coordinates ≈ T[
  1  -0   0  -0;
  0   0   1  -0;
 -3   3  -2  -1;
  2  -2   1   1;
] # from https://defelement.org/elements/examples/interval-hermite-3.html

@test get_face_own_dofs(hermite1) == [ [1, 3],    [2, 4],    [] ]


hermite2 = HermiteRefFE(T, TRI)
dofs = get_dof_basis(hermite2)
shapefuns = get_shapefuns(hermite2)
@test norm(evaluate(dofs, shapefuns) - I) <= 1.e-15

P3_2D_monoms = FEEC_poly_basis(Val(2), T, 3, 0, :P, Monomial)
shapefuns_monom_coordinates = inv(evaluate(dofs, P3_2D_monoms))
@test shapefuns_monom_coordinates ≈ T[
   1  -0   0    0   0   0   0  -0   0  -0;
   0   0   0    0   1   0   0  -0   0  -0;
  -3   3   0    0  -2   0  -1  -0   0  -0;
   2  -2   0    0   1   0   1  -0   0  -0;
   0   0   0    0   0   1   0  -0   0  -0;
 -13  -7  -7   27  -3  -3   2  -1  -1   2;
  13   7   7  -27   3   2  -2   2   1  -2;
  -3   0   3   -0   0  -2   0   0   0  -1;
  13   7   7  -27   2   3  -2   1   2  -2;
   2   0  -2    0   0   1   0   0   0   1;
] # from https://defelement.org/elements/examples/triangle-hermite-3.html

@test get_face_own_dofs(hermite2) == [
#  φ  ∂x(φ) ∂y(φ)  at v1
  [1, 5,    6],
#  φ  ∂x(φ) ∂y(φ)  at v2
  [2, 7,    8],
#  φ  ∂x(φ) ∂y(φ)  at v3
  [3, 9,    10],
  [], [], [],
#  φ at TRI barycenter
  [4]
]
@test get_face_dofs(hermite2) == [
  [1, 5, 6], [2, 7, 8], [3, 9, 10],
  [1, 5, 6, 2, 7, 8], [1, 5, 6, 3, 9, 10], [2, 7, 8, 3, 9, 10],
  [1, 5, 6, 2, 7, 8, 3, 9, 10, 4]
]

V = SymTracelessTensorValue{2,T}
hermite1_q = HermiteRefFE(V, SEGMENT)
dofs = get_dof_basis(hermite1_q)
shapefuns = get_shapefuns(hermite1_q)
@test norm(evaluate(dofs, shapefuns) - I) <= 1.e-15

P3_1D_monoms = FEEC_poly_basis(Val(1), V, 3, 0, :P, Monomial; cart_prod=true)
shapefuns_monom_coordinates = inv(evaluate(dofs, P3_1D_monoms))
shapefuns_monom_coordinates ≈ T[
# v1  v2  v1  v2  v1   v1   v2   v2
# φ₁  φ₁  φ₂  φ₂ ∂xφ₁ ∂xφ₂ ∂xφ₁ ∂xφ₂
  1   0  -0  -0   0    0    0   -0;
  0   0   1  -0   0    0    0   -0;
  0   0   0  -0   1    0    0   -0;
  0   0   0   0   0    1    0   -0;
 -3   3   0   0  -2   -0   -1   -0;
  0   0  -3   3   0   -2    0   -1;
  2  -2   0   0   1    0    1   -0;
  0   0   2  -2   0    1    0    1;
]

@test get_face_own_dofs(hermite1_q) == Vector{Int}[
#  φ₁ φ₂ ∂x(φ₁) ∂x(φ₂)  at v1
  [1, 3,  5,    6],
#  φ₁ φ₂ ∂x(φ₁) ∂x(φ₂)  at v2
  [2, 4,  7,    8],
  []
]


V = VectorValue{2,T}
hermite3_v = HermiteRefFE(V, TET)
dofs = get_dof_basis(hermite3_v)
shapefuns = get_shapefuns(hermite3_v)
@test norm(evaluate(dofs, shapefuns) - I) <= 1.e-15
@test get_face_own_dofs(hermite3_v) == Vector{Int}[
#  φ₁ φ₂ ∂x(φ₁) ∂x(φ₂)  at v3
  [1, 9,  17, 18, 19, 20, 21, 22],
  [2, 10, 23, 24, 25, 26, 27, 28],
  [3, 11, 29, 30, 31, 32, 33, 34],
  [4, 12, 35, 36, 37, 38, 39, 40],
  [], [], [], [], [], [],
#  φ₁ φ₂  at face 1 barycenter
  [5, 13],
#  φ₁ φ₂  at face 2 barycenter
  [6, 14],
#  φ₁ φ₂  at face 3 barycenter
  [7, 15],
#  φ₁ φ₂  at face 3 barycenter
  [8, 16],
  []
]

# Physical mapping of the element TODO

end
