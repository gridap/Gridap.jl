module DofsTests

using Test
using Gridap
using Gridap.ReferenceFEs
using Gridap.Fields
using FillArrays

T = Float64

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

end # module
