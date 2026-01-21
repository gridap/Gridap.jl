module DofsTests

using Test
using Gridap
using Gridap.ReferenceFEs
using Gridap.Fields
using FillArrays

T = Float64

# Linear combination of Dof bases

reffe0 = ReferenceFE(SEGMENT, Lagrangian(), Float64, 2)
reffe1 = ReferenceFE(SEGMENT, Lagrangian(), Float64, 3)

dofs  = get_dof_basis(reffe0)                  # length 3
@test dofs    isa AbstractVector{<:Dof}

dofs_lc = linear_combination(Eye(3,6), dofs)   # length 6
@test dofs_lc isa AbstractVector{<:Dof}
@test length(dofs_lc) == 6

basis = get_prebasis(reffe1)                   # length 4
@test basis    isa AbstractVector{<:Field}
basis_lc = linear_combination(Eye(4,5), basis) # length 5

@test basis_lc isa AbstractVector{<:Field}
@test length(basis_lc) == 5

M34 = evaluate(dofs,    basis   )
M35 = evaluate(dofs,    basis_lc)
M64 = evaluate(dofs_lc, basis   )
M65 = evaluate(dofs_lc, basis_lc)

@test size(M34) == (3, 4)
@test size(M35) == (3, 5)
@test size(M64) == (6, 4)
@test size(M65) == (6, 5)
@test M35[1:3,1:4] == M34
@test M64[1:3,1:4] == M34
@test M65[1:3,1:4] == M34

# Concatenation of Dof bases

P = Point{2,Float64}
V = VectorValue{3,Float64}
φ = ConstantField(V(1,2,3))
φ2= ConstantField(V(4,5,6))
φφ2 = [φ, φ2]

dofs0 = LagrangianDofBasis(V, P[rand(P) for _ in 1:0])
dofs_0 = vcat(dofs0)
@test length(dofs_0) == length(dofs0)
@test evaluate(dofs_0, φ) == evaluate(dofs0, φ)

dofs1 = LagrangianDofBasis(V, [rand(P) for _ in 1:2])
dofs_1 = vcat(dofs1)
@test length(dofs_1) == length(dofs1)
@test evaluate(dofs_1, φ) == evaluate(dofs1, φ)
@test evaluate(dofs_1, φφ2) == evaluate(dofs1, φφ2)

dofs_01 = vcat(dofs0, dofs1)
dofs_10 = vcat(dofs1, dofs0)
@test length(dofs1) == length(dofs_01) == length(dofs_10)
@test evaluate(dofs1, φ) == evaluate(dofs_01, φ) == evaluate(dofs_10, φ)
@test evaluate(dofs1, φφ2) == evaluate(dofs_01, φφ2) == evaluate(dofs_10, φφ2)

dofs2 = get_dof_basis(CrouzeixRaviartRefFE(V,TRI,1))
dofs_12 = vcat(dofs1, dofs2)
@test length(dofs_12) == length(dofs1) + length(dofs2)
dofs1_at_φ = evaluate(dofs1, φ)
dofs2_at_φ = evaluate(dofs2, φ)
@test evaluate(dofs_12, φ) == vcat(dofs1_at_φ, dofs2_at_φ)
dofs1_at_φφ2 = evaluate(dofs1, φφ2)
dofs2_at_φφ2 = evaluate(dofs2, φφ2)
@test evaluate(dofs_12, φφ2) == vcat(dofs1_at_φφ2, dofs2_at_φφ2)

dofs_00200100 = vcat(dofs0, dofs0, dofs2, dofs0, dofs0, dofs1, dofs0, dofs0)
@test evaluate(dofs_00200100, φ)   == vcat(dofs2_at_φ, dofs1_at_φ)
@test evaluate(dofs_00200100, φφ2) == vcat(dofs2_at_φφ2, dofs1_at_φφ2)

end # module
