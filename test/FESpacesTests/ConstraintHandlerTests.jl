using Test
using LinearAlgebra
using SparseArrays
using Gridap
using Gridap.FESpaces
using Gridap.Arrays: Table, datarange

# -----------------------------------------------------------------------
# Construction and basic queries
# -----------------------------------------------------------------------

ch = ConstraintHandler(10)
@test num_constrained_dofs(ch) == 0
@test num_free_dofs(ch) == 10
@test !is_constrained(ch, 3)

add_constraint!(ch, 3, [1, 2], [0.5, 0.5])
@test is_constrained(ch, 3)
@test num_constrained_dofs(ch) == 1
@test num_free_dofs(ch) == 9

# -----------------------------------------------------------------------
# Dirichlet (pure offset)
# -----------------------------------------------------------------------

ch = ConstraintHandler(5)
add_constraint!(ch, 2, 3.0)
close!(ch)
x = zeros(5)
apply_constraints!(x, ch)
@test x[2] ≈ 3.0

# -----------------------------------------------------------------------
# Homogeneous linear constraint
# -----------------------------------------------------------------------

ch = ConstraintHandler(5)
add_constraint!(ch, 3, [1, 2], [0.5, 0.5])
close!(ch)
x = [2.0, 4.0, 0.0, 0.0, 0.0]
apply_constraints!(x, ch)
@test x[3] ≈ 3.0

# -----------------------------------------------------------------------
# Affine constraint (linear + offset)
# -----------------------------------------------------------------------

ch = ConstraintHandler(4)
add_constraint!(ch, 4, [1], [2.0], 1.0)
close!(ch)
x = [3.0, 0.0, 0.0, 0.0]
apply_constraints!(x, ch)
@test x[4] ≈ 7.0

# -----------------------------------------------------------------------
# set_offset!
# -----------------------------------------------------------------------

ch = ConstraintHandler(3)
add_constraint!(ch, 2, 1.0)
set_offset!(ch, 2, 5.0)
close!(ch)
x = zeros(3)
apply_constraints!(x, ch)
@test x[2] ≈ 5.0

# -----------------------------------------------------------------------
# Chain resolution: one level
# -----------------------------------------------------------------------

ch = ConstraintHandler(5)
add_constraint!(ch, 2, [1], [1.0])
add_constraint!(ch, 3, [2], [1.0])  # x3 = x2 = x1
close!(ch)

line3 = ch.constraints[findfirst(l -> l.dof == 3, ch.constraints)]
@test length(line3.masters) == 1
@test line3.masters[1] == 1
@test line3.coeffs[1] ≈ 1.0

x = [7.0, 0.0, 0.0, 0.0, 0.0]
apply_constraints!(x, ch)
@test x[2] ≈ 7.0
@test x[3] ≈ 7.0

# -----------------------------------------------------------------------
# Chain resolution: weighted chain
# -----------------------------------------------------------------------

ch = ConstraintHandler(4)
add_constraint!(ch, 2, [1], [2.0])
add_constraint!(ch, 3, [2], [3.0])  # x3 = 6*x1
close!(ch)

line3 = ch.constraints[findfirst(l -> l.dof == 3, ch.constraints)]
@test length(line3.masters) == 1
@test line3.masters[1] == 1
@test line3.coeffs[1] ≈ 6.0

x = [1.0, 0.0, 0.0, 0.0]
apply_constraints!(x, ch)
@test x[2] ≈ 2.0
@test x[3] ≈ 6.0

# -----------------------------------------------------------------------
# Chain resolution: offset propagation
# -----------------------------------------------------------------------

# x2 = x1 + 1,  x3 = x2 + 2  →  x3 = x1 + 3
ch = ConstraintHandler(5)
add_constraint!(ch, 2, [1], [1.0], 1.0)
add_constraint!(ch, 3, [2], [1.0], 2.0)
close!(ch)

x = [10.0, 0.0, 0.0, 0.0, 0.0]
apply_constraints!(x, ch)
@test x[2] ≈ 11.0
@test x[3] ≈ 13.0

# -----------------------------------------------------------------------
# Duplicate master entry merging
# -----------------------------------------------------------------------

ch = ConstraintHandler(3)
add_constraint!(ch, 3, [1, 1, 2], [0.3, 0.2, 0.5])
close!(ch)

masters = Dict(zip(ch.constraints[1].masters, ch.constraints[1].coeffs))
@test masters[1] ≈ 0.5
@test masters[2] ≈ 0.5

# -----------------------------------------------------------------------
# Near-zero coefficient stripping
# -----------------------------------------------------------------------

ch = ConstraintHandler(3)
add_constraint!(ch, 3, [1, 2], [1.0, 1e-15])
close!(ch)
@test length(ch.constraints[1].masters) == 1
@test ch.constraints[1].masters[1] == 1

# -----------------------------------------------------------------------
# Idempotent close!
# -----------------------------------------------------------------------

ch = ConstraintHandler(3)
add_constraint!(ch, 2, [1], [1.0])
close!(ch)
close!(ch)
@test ch.closed
@test num_constrained_dofs(ch) == 1

# -----------------------------------------------------------------------
# Cycle detection
# -----------------------------------------------------------------------

ch = ConstraintHandler(5)
add_constraint!(ch, 2, [3], [1.0])
add_constraint!(ch, 3, [2], [1.0])  # cycle: 2 ↔ 3
@test_throws ErrorException close!(ch)
@test occursin("cycle", sprint(showerror, try close!(ch) catch e e end))

# -----------------------------------------------------------------------
# merge!: no conflict
# -----------------------------------------------------------------------

ch1 = ConstraintHandler(6)
add_constraint!(ch1, 1, [3], [1.0])
ch2 = ConstraintHandler(6)
add_constraint!(ch2, 2, [4, 5], [0.5, 0.5])
merge!(ch1, ch2)
close!(ch1)
@test is_constrained(ch1, 1)
@test is_constrained(ch1, 2)
x = [0.0, 0.0, 7.0, 2.0, 4.0, 0.0]
apply_constraints!(x, ch1)
@test x[1] ≈ 7.0
@test x[2] ≈ 3.0

# -----------------------------------------------------------------------
# merge!: chain across two handlers
# -----------------------------------------------------------------------

ch1 = ConstraintHandler(4)
add_constraint!(ch1, 3, [2], [1.0])  # x3 = x2
ch2 = ConstraintHandler(4)
add_constraint!(ch2, 2, [1], [2.0])  # x2 = 2*x1
merge!(ch1, ch2)
close!(ch1)  # should resolve x3 = x2 = 2*x1

line3 = ch1.constraints[findfirst(l -> l.dof == 3, ch1.constraints)]
@test length(line3.masters) == 1
@test line3.masters[1] == 1
@test line3.coeffs[1] ≈ 2.0

# -----------------------------------------------------------------------
# merge!: conflict raises error by default
# -----------------------------------------------------------------------

ch1 = ConstraintHandler(4)
add_constraint!(ch1, 2, [1], [1.0])
ch2 = ConstraintHandler(4)
add_constraint!(ch2, 2, [3], [1.0])  # conflict on dof 2
@test_throws ErrorException merge!(ch1, ch2)

# -----------------------------------------------------------------------
# merge!: user-supplied conflict resolver
# -----------------------------------------------------------------------

ch1 = ConstraintHandler(4)
add_constraint!(ch1, 2, [1], [1.0])
ch2 = ConstraintHandler(4)
add_constraint!(ch2, 2, [3], [1.0])

# Average the two constraints: x2 = 0.5*x1 + 0.5*x3
resolve = (existing, incoming) ->
  ConstraintLine(existing.dof,
    vcat(existing.masters, incoming.masters),
    vcat(existing.coeffs, incoming.coeffs) .* 0.5,
    0.0)

merge!(ch1, ch2; on_conflict=resolve)
close!(ch1)
x = [4.0, 0.0, 6.0, 0.0]
apply_constraints!(x, ch1)
@test x[2] ≈ 5.0

# -----------------------------------------------------------------------
# apply_constraints (system reduction): pure Dirichlet
# -----------------------------------------------------------------------

# 2 dofs: x1 free, x2 = 3.0 (Dirichlet)
# Ã = [2], b̃ = [5 - 1*3] = [2]  →  u = [1]  →  x = [1, 3]
A = sparse([2.0 1.0; 1.0 2.0])
b = [5.0, 6.0]
ch = ConstraintHandler(2)
add_constraint!(ch, 2, 3.0)
close!(ch)

Ã, b̃ = apply_constraints(A, b, ch)
C, g  = get_matrix(ch)

@test size(Ã) == (1, 1)
u = Ã \ b̃
x = C * u + g
@test x[1] ≈ 1.0
@test x[2] ≈ 3.0

# -----------------------------------------------------------------------
# apply_constraints (system reduction): linear constraint
# -----------------------------------------------------------------------

# 3 dofs: x3 = 0.5*x1 + 0.5*x2, A = I
# b = A * [1, 2, 1.5] = [1, 2, 1.5] (consistent with constraint)
n = 3
A = sparse(1.0I, n, n)
x_target = [1.0, 2.0, 1.5]
b = A * x_target

ch = ConstraintHandler(n)
add_constraint!(ch, 3, [1, 2], [0.5, 0.5])
close!(ch)

Ã, b̃ = apply_constraints(A, b, ch)
C, g  = get_matrix(ch)

@test size(Ã) == (2, 2)
u = Ã \ b̃
x = C * u + g
@test x[1] ≈ 1.0
@test x[2] ≈ 2.0
@test x[3] ≈ 1.5

# -----------------------------------------------------------------------
# apply_constraints (system reduction): full system round-trip
# -----------------------------------------------------------------------

# Poisson-like system with a known analytical solution.
# 4 dofs: symmetric tridiagonal stiffness with constraints:
#   x1 = 0 (Dirichlet left)
#   x4 = 1 (Dirichlet right)
n = 4
A = sparse([2.0 -1.0  0.0  0.0;
           -1.0  2.0 -1.0  0.0;
            0.0 -1.0  2.0 -1.0;
            0.0  0.0 -1.0  2.0])
b = zeros(n)

ch = ConstraintHandler(n)
add_constraint!(ch, 1, 0.0)
add_constraint!(ch, 4, 1.0)
close!(ch)

Ã, b̃ = apply_constraints(A, b, ch)
C, g  = get_matrix(ch)

@test size(Ã) == (2, 2)
u = Ã \ b̃
x = C * u + g

# Linear interpolation between 0 and 1 is exact for this 1D Poisson.
@test x[1] ≈ 0.0   atol=1e-12
@test x[2] ≈ 1/3   atol=1e-12
@test x[3] ≈ 2/3   atol=1e-12
@test x[4] ≈ 1.0   atol=1e-12

# -----------------------------------------------------------------------
# Error cases
# -----------------------------------------------------------------------

ch = ConstraintHandler(5)
close!(ch)
@test_throws AssertionError add_constraint!(ch, 1, 0.0)

ch = ConstraintHandler(5)
add_constraint!(ch, 2, [1], [1.0])
@test_throws ErrorException add_constraint!(ch, 2, [3], [1.0])

ch = ConstraintHandler(5)
@test_throws AssertionError add_constraint!(ch, 0, 0.0)
@test_throws AssertionError add_constraint!(ch, 6, 0.0)

ch = ConstraintHandler(5)
@test_throws AssertionError add_constraint!(ch, 3, [3], [1.0])

ch = ConstraintHandler(3)
add_constraint!(ch, 2, [1], [1.0])
@test_throws AssertionError apply_constraints!(zeros(3), ch)

ch1 = ConstraintHandler(3)
close!(ch1)
ch2 = ConstraintHandler(3)
@test_throws AssertionError merge!(ch1, ch2)

ch1 = ConstraintHandler(4)
ch2 = ConstraintHandler(5)
@test_throws AssertionError merge!(ch1, ch2)

# -----------------------------------------------------------------------
# constraint_tables — requires a real FESpace
# -----------------------------------------------------------------------

# Space: 6 free + 3 Dirichlet dofs (2×2 Cartesian mesh, Lagrangian order 1,
# Dirichlet on left edge x=0).
# ConstraintHandler DoF convention: 1..6 free, 7..9 Dirichlet.
# Gridap signed-dof convention:     1..6 free, -1..-3 Dirichlet.

const _model  = CartesianDiscreteModel((0,1,0,1), (2,2))
const _reffe  = ReferenceFE(lagrangian, Float64, 1)
let labels = get_face_labeling(_model)
  add_tag_from_tags!(labels, "left", ["tag_7", "tag_1", "tag_3"])
end
const _space  = FESpace(_model, _reffe; dirichlet_tags="left")

# Free master
ch = ConstraintHandler(num_free_dofs(_space) + num_dirichlet_dofs(_space))
add_constraint!(ch, 3, [1, 2], [0.5, 0.5])
close!(ch)
sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs = constraint_tables(_space, ch)
@test sDOF_to_dof == [3]
@test collect(sDOF_to_dofs[1])   == [1, 2]
@test collect(sDOF_to_coeffs[1]) ≈  [0.5, 0.5]

# Dirichlet master (ConstraintHandler DoF 7 → signed dof -1)
ch = ConstraintHandler(num_free_dofs(_space) + num_dirichlet_dofs(_space))
add_constraint!(ch, 3, [7], [1.0])
close!(ch)
sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs = constraint_tables(_space, ch)
@test sDOF_to_dof == [3]
@test collect(sDOF_to_dofs[1]) == [-1]

# Slave is a Dirichlet dof (ConstraintHandler DoF 7 → signed dof -1)
ch = ConstraintHandler(num_free_dofs(_space) + num_dirichlet_dofs(_space))
add_constraint!(ch, 7, [1], [1.0])
close!(ch)
sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs = constraint_tables(_space, ch)
@test sDOF_to_dof == [-1]
@test collect(sDOF_to_dofs[1]) == [1]

# Chain resolved before output: x3 = x2 = x1
ch = ConstraintHandler(num_free_dofs(_space) + num_dirichlet_dofs(_space))
add_constraint!(ch, 2, [1], [1.0])
add_constraint!(ch, 3, [2], [1.0])
close!(ch)
sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs = constraint_tables(_space, ch)
i3 = findfirst(==(3), sDOF_to_dof)
@test collect(sDOF_to_dofs[i3])   == [1]
@test collect(sDOF_to_coeffs[i3]) ≈  [1.0]

# Round-trip: constraint_tables → ConstraintHandler
ch_orig = ConstraintHandler(num_free_dofs(_space) + num_dirichlet_dofs(_space))
add_constraint!(ch_orig, 3, [1, 2], [0.5, 0.5])
add_constraint!(ch_orig, 5, [4],    [1.0])
close!(ch_orig)
sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs = constraint_tables(_space, ch_orig)
ch_rt = ConstraintHandler(sDOF_to_dof, sDOF_to_dofs, sDOF_to_coeffs, _space)
close!(ch_rt)
@test num_constrained_dofs(ch_rt) == num_constrained_dofs(ch_orig)
@test is_constrained(ch_rt, 3)
@test is_constrained(ch_rt, 5)
