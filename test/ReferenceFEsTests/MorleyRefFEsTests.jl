module MorleyRefFEsTest

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

using Gridap.Geometry: Point
using LinearAlgebra

T = Float64

# Morley DOF tests TODO


# RefFE tests

morley = MorleyRefFE(T, TRI, 2)
dofs = get_dof_basis(morley)
shapefuns = get_shapefuns(morley)
@test norm(evaluate(dofs, shapefuns) - I) <= 1.e-14

o_sq2 = 0.7071067811865475 # 1/sqrt2
@test evaluate(dofs, MonomialBasis{2}(T,2))[4:6,:] == T[
  0.0 0.0 0.0 -1.0 -0.5 -0.25 0.0 0.0 0.0;
  0.0 -1.0 0.0 0.0 -0.5 0.0 0.0 -0.25 0.0;
  0.0 o_sq2 o_sq2 o_sq2 o_sq2 3*o_sq2/4 o_sq2 3*o_sq2/4 o_sq2/2
]

#c = return_cache(dofs, shapefuns)
#@btime evaluate!($c, $dofs, $shapefuns)

#@profview for _ in 1:500000 evaluate!(c, dofs, shapefuns) end
#@profview_allocs for _ in 1:500000 evaluate!(c, dofs, shapefuns) end

@test get_face_own_dofs(morley) == Array{Int}[
#  φ at vertices v1, v2, v3
  [1], [2], [3],
#  ∂ₙφ deᵢ at edges e₁,e₂,e₂ midpoints
  [4], [5], [6],
  []
]
@test get_face_dofs(morley) == Array{Int}[
  [1], [2], [3],
  [1, 2, 4],
  [1, 3, 5],
  [2, 3, 6],
  [1, 2, 3, 4, 5, 6]
]

V = VectorValue{2,T}
morley_q = MorleyRefFE(V, TRI, 2)
dofs = get_dof_basis(morley_q)
shapefuns = get_shapefuns(morley_q)
@test norm(evaluate(dofs, shapefuns) - I) <= 1.e-13

o_sq2 = 0.7071067811865475 # 1/sqrt2
@test evaluate(dofs, MonomialBasis{2}(V,1))[7:12,:] == T[
  0.0 0.0  0.0   0.0 -1.0    0.0   -0.5   0.0 ;
  0.0 0.0  0.0   0.0  0.0   -1.0    0.0  -0.5 ;
  0.0 0.0 -1.0   0.0  0.0    0.0   -0.5   0.0 ;
  0.0 0.0  0.0  -1.0  0.0    0.0    0.0  -0.5 ;
  0.0 0.0 o_sq2  0.0  o_sq2   0.0  o_sq2  0.0
  0.0 0.0  0.0  o_sq2  0.0   o_sq2  0.0  o_sq2
]


@test get_face_own_dofs(morley_q) == Array{Int}[
#  φ₁ φ₂  at vertices v1, v2, v3
  [1, 2], [3, 4],  [5, 6],
#  ∂ₙφᵢ  at edges e₁,e₂,e₂ midpoints
  [7, 8], [9, 10], [11, 12],
  []
]


# Physical mapping of the element TODO

# Reference and physical polytope
K̂ = GeneralPolytope{2}(TRI);

# Φ: K̂ -> K
Φ = affine_map(TensorValue(1., .2, .3, 1.), Point(1,2.)); # Translated, JΦ = I₂
# It's inverse: should'nt be used in practice because the only physical points
# at which the cell_field are usually evaluated are the mapped reference point
Φ_inv = inverse_map(Φ);

v = evaluate(Φ, get_vertex_coordinates(TRI))
K = GeneralPolytope{2}(TRI, v);

V = T # VectorValue{2,T}
ref_fe = MorleyRefFE(V, TRI, 2; vertices=get_vertex_coordinates(K̂));
phy_fe = MorleyRefFE(V, TRI, 2; vertices=get_vertex_coordinates(K));

# Reference and physical shape function, computed explicitly
σ̂ = get_dof_basis(ref_fe); φ̂ = get_shapefuns(ref_fe); norm(evaluate(σ̂,φ̂)-I)
σ = get_dof_basis(phy_fe); φ = get_shapefuns(phy_fe); norm(evaluate(σ,φ)-I)

# Let's try to construct φ and σ from φ̂, σ̂, Φ, and σ to find the COB.
C_σ_Φφ̂ = evaluate(σ, φ̂ .∘ Φ_inv); # manual shapefun pushforward
round.(inv(C_σ_Φφ̂); digits=2) # DoF change matrix
Jt = Array(Φ.gradient)
round.(C_σ_Φφ̂'; digits=5) # shapefuns change matrix
invJt = Array(Φ_inv.gradient')

φ_pushed  = linear_combination(inv(C_σ_Φφ̂), φ̂ .∘ Φ_inv);
norm(evaluate(σ, φ_pushed)-I)
σ_pushed = linear_combination(transpose(C_σ_Φφ̂), σ̂);

norm(evaluate(σ_pushed,  φ_pushed.∘Φ)-I) # manual dof pullback

end
