module ConstraintHandlerTests

using Test
using LinearAlgebra
using SparseArrays
using Gridap
import Gridap.Geometry
using Gridap.FESpaces
using Gridap.FESpaces: select_constraint_slaves
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

_model  = CartesianDiscreteModel((0,1,0,1), (2,2))
_reffe  = ReferenceFE(lagrangian, Float64, 1)
let labels = get_face_labeling(_model)
  add_tag_from_tags!(labels, "left", ["tag_7", "tag_1", "tag_3"])
end
_space  = FESpace(_model, _reffe; dirichlet_tags="left")

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


# -----------------------------------------------------------------------
# merge_slave_constraint_tables / close_slave_constraint_tables
# -----------------------------------------------------------------------

using Gridap.FESpaces: merge_slave_constraint_tables, close_slave_constraint_tables

# Helper: slave_dof → Dict(master_dof → coeff), order-independent comparison.
_cdict(dofs, dtable, ctable) = Dict(
  dofs[i] => Dict(zip(collect(dtable[i]), collect(ctable[i])))
  for i in eachindex(dofs))

_n_dofs = num_free_dofs(_space) + num_dirichlet_dofs(_space)

# Merge: disjoint sets → fast path, all slaves present in output.
let
  ch1 = ConstraintHandler(_n_dofs); add_constraint!(ch1, 3, [1,2], [0.5,0.5]); close!(ch1)
  ch2 = ConstraintHandler(_n_dofs); add_constraint!(ch2, 5, [4],   [1.0]);     close!(ch2)
  d1,t1,c1,o1 = constraint_tables(_space, ch1)
  d2,t2,c2,o2 = constraint_tables(_space, ch2)
  out_dof, out_dofs, out_coeffs, _ = merge_slave_constraint_tables(
    _space, d1,t1,c1, d2,t2,c2, o1,o2)
  @test Set(out_dof) == Set([3, 5])
  @test _cdict(out_dof, out_dofs, out_coeffs)[3] == Dict(1 => 0.5, 2 => 0.5)
  @test _cdict(out_dof, out_dofs, out_coeffs)[5] == Dict(4 => 1.0)
end

# Merge: conflict → on_conflict (first-wins keeps s1 row).
let
  ch1 = ConstraintHandler(_n_dofs); add_constraint!(ch1, 3, [1], [1.0]); close!(ch1)
  ch2 = ConstraintHandler(_n_dofs); add_constraint!(ch2, 3, [2], [1.0]); close!(ch2)
  d1,t1,c1,o1 = constraint_tables(_space, ch1)
  d2,t2,c2,o2 = constraint_tables(_space, ch2)
  first_wins = (dof,md1,mc1,mo1,md2,mc2,mo2) -> (collect(md1), collect(mc1), mo1)
  out_dof, out_dofs, out_coeffs, _ = merge_slave_constraint_tables(
    _space, d1,t1,c1, d2,t2,c2, o1,o2; on_conflict=first_wins)
  @test _cdict(out_dof, out_dofs, out_coeffs)[3] == Dict(1 => 1.0)
end

# Hard test: path A (merge ConstraintHandlers → close → tables) must equal
# path B (tables per handler → merge_slave → close_slave).
#
# ch1: dof3=dof2, dof5=0.5*dof1+0.5*dof4
# ch2: dof2=dof1, dof5=dof4  (conflict on dof5; ch1/s1 wins)
#
# After merge+close: dof2=dof1, dof3=dof1 (chain dof3→dof2→dof1), dof5=0.5*dof1+0.5*dof4.
# Path B exercises close_slave_constraint_tables on the cross-handler chain dof3→dof2.
let
  first_wins_ch  = (existing, _) -> existing
  first_wins_tbl = (dof,md1,mc1,mo1,md2,mc2,mo2) -> (collect(md1), collect(mc1), mo1)

  # Path A
  ch_a1 = ConstraintHandler(_n_dofs)
  add_constraint!(ch_a1, 3, [2],   [1.0])
  add_constraint!(ch_a1, 5, [1,4], [0.5, 0.5])
  ch_a2 = ConstraintHandler(_n_dofs)
  add_constraint!(ch_a2, 2, [1], [1.0])
  add_constraint!(ch_a2, 5, [4], [1.0])
  merge!(ch_a1, ch_a2; on_conflict=first_wins_ch)
  close!(ch_a1)
  dof_A, dofs_A, coeffs_A, _ = constraint_tables(_space, ch_a1)

  # Path B
  ch_b1 = ConstraintHandler(_n_dofs)
  add_constraint!(ch_b1, 3, [2],   [1.0])
  add_constraint!(ch_b1, 5, [1,4], [0.5, 0.5])
  close!(ch_b1)
  d1,t1,c1,o1 = constraint_tables(_space, ch_b1)

  ch_b2 = ConstraintHandler(_n_dofs)
  add_constraint!(ch_b2, 2, [1], [1.0])
  add_constraint!(ch_b2, 5, [4], [1.0])
  close!(ch_b2)
  d2,t2,c2,o2 = constraint_tables(_space, ch_b2)

  dof_B, dofs_B, coeffs_B, offs_B = merge_slave_constraint_tables(
    _space, d1,t1,c1, d2,t2,c2, o1,o2; on_conflict=first_wins_tbl)
  dof_B, dofs_B, coeffs_B, _ = close_slave_constraint_tables(
    _space, dof_B, dofs_B, coeffs_B, offs_B)

  @test _cdict(dof_A, dofs_A, coeffs_A) == _cdict(dof_B, dofs_B, coeffs_B)
end

# ===========================================================================
# ConstraintHandler — cell-array constructor
# ===========================================================================

# Single constraint per cell (vector form)
let
  n_dofs = 8
  cell_slave_ids  = [3, 5, 7]
  cell_master_ids = [[1, 2], [4, 6], [2, 4, 6]]
  cell_coeffs     = [[0.5, 0.5], [1.0, 0.0], [0.0, 0.5, 0.5]]
  cell_offsets    = [0.0, 1.0, -0.5]

  ch = ConstraintHandler(n_dofs, cell_slave_ids, cell_master_ids, cell_coeffs, cell_offsets)

  @test num_constrained_dofs(ch) == 3
  for dof in [3, 5, 7]; @test is_constrained(ch, dof); end
  line5 = ch.constraints[ch.dof_to_constraint[5]]
  @test line5.coeffs ≈ [1.0, 0.0]
  @test line5.offset ≈  1.0
end

# Multiple constraints per cell (matrix form)
let
  n_dofs = 6
  cell_slave_ids2  = [[3, 4]]
  cell_master_ids2 = [[1, 2]]
  cell_coeffs2     = [[ 0.5  0.5; 0.0  1.0]]
  cell_offsets2    = [[1.0, -1.0]]

  ch2 = ConstraintHandler(n_dofs, cell_slave_ids2, cell_master_ids2, cell_coeffs2, cell_offsets2)
  @test num_constrained_dofs(ch2) == 2
  @test ch2.constraints[ch2.dof_to_constraint[3]].coeffs ≈ [0.5, 0.5]
  @test ch2.constraints[ch2.dof_to_constraint[4]].offset ≈ -1.0
  close!(ch2)
  @test ch2.closed
end

# on_conflict: error by default
@test_throws ErrorException ConstraintHandler(4, [2, 2], [[1], [3]], [[1.0], [2.0]], [0.0, 0.0])

# on_conflict: user resolver (last wins)
let
  ch_last = ConstraintHandler(4, [2, 2], [[1], [3]], [[1.0], [2.0]], [0.0, 0.0];
                               on_conflict = (e, i) -> i)
  @test ch_last.constraints[ch_last.dof_to_constraint[2]].masters == [3]
end

# ===========================================================================
# select_constraint_slaves
# ===========================================================================

# Basic: 3 constraints, trivial selection
let
  dof_ids = [[1, 2], [3], [4, 5]]
  coeffs  = [[1.0, 2.0], [3.0], [1.0, 1.0]]
  offsets = [4.0, 6.0, 0.0]

  slave_ids, masters, coeffs_t, offs = select_constraint_slaves(5, dof_ids, coeffs, offsets)

  @test slave_ids[1] == 2
  @test collect(masters[1]) == [1]
  @test collect(coeffs_t[1]) ≈ [-0.5]
  @test offs[1] ≈ 2.0

  @test slave_ids[2] == 3
  @test collect(masters[2]) == []
  @test offs[2] ≈ 2.0

  s3 = slave_ids[3]
  m3 = s3 == 4 ? 5 : 4
  @test s3 ∈ (4, 5)
  @test collect(masters[3]) == [m3]
  @test collect(coeffs_t[3]) ≈ [-1.0]
  @test offs[3] ≈ 0.0
end

# Shared DOF: distinct slaves guaranteed
let
  dof_ids = [[1, 2], [2, 3]]
  coeffs  = [[1.0, 1.0], [1.0, 1.0]]
  offsets = [0.0, 0.0]

  slave_ids, masters, coeffs_t, offs = select_constraint_slaves(3, dof_ids, coeffs, offsets)

  @test length(slave_ids) == 2
  @test slave_ids[1] != slave_ids[2]
  for i in 1:2
    @test collect(coeffs_t[i]) ≈ [-1.0]
    @test offs[i] ≈ 0.0
  end
end

# Connectivity tiebreaker
let
  dof_ids = [[1, 2], [2, 3], [2, 4], [2, 5]]
  coeffs  = [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0], [1.0, 1.0]]
  offsets = [0.0, 0.0, 0.0, 0.0]

  slave_ids, _, _, _ = select_constraint_slaves(5, dof_ids, coeffs, offsets)
  @test slave_ids[1] == 1
  @test 2 ∈ slave_ids
end

# Zero coefficients are ignored
let
  dof_ids = [[1, 2, 3]]
  coeffs  = [[1.0, 0.0, 2.0]]
  offsets = [6.0]

  slave_ids, masters, coeffs_t, offs = select_constraint_slaves(3, dof_ids, coeffs, offsets)

  @test slave_ids[1] == 3
  @test collect(masters[1]) == [1]
  @test collect(coeffs_t[1]) ≈ [-0.5]
  @test offs[1] ≈ 3.0
end

# Infeasible: no slave candidate → error
let
  dof_ids = [[1], [2], [1, 2]]
  coeffs  = [[1.0], [1.0], [1.0, 1.0]]
  offsets = [0.0, 0.0, 0.0]
  @test_throws Exception select_constraint_slaves(2, dof_ids, coeffs, offsets)
end

# Multiple constraints per cell (matrix form)
let
  cell_dof_ids = [[1, 2, 3]]
  cell_coeffs  = [[1.0  1.0  0.0;
                   0.0  1.0  2.0]]
  cell_offsets = [[0.0, 6.0]]

  slave_ids, masters, coeffs_t, offs =
    select_constraint_slaves(3, cell_dof_ids, cell_coeffs, cell_offsets)

  @test length(slave_ids) == 2
  @test slave_ids[1] != slave_ids[2]
  @test slave_ids[1] == 1
  @test collect(masters[1]) == [2]
  @test collect(coeffs_t[1]) ≈ [-1.0]
  @test offs[1] ≈ 0.0
  @test slave_ids[2] == 3
  @test collect(masters[2]) == [2]
  @test collect(coeffs_t[2]) ≈ [-0.5]
  @test offs[2] ≈ 3.0
end

# ===========================================================================
# Constructor 1 — cell, bilinear: DG P2 \ P1
# ===========================================================================

let
  model = CartesianDiscreteModel((0,1), (1,))
  V = FESpace(model, ReferenceFE(lagrangian, Float64, 2); conformity=:L2)
  W = FESpace(model, ReferenceFE(lagrangian, Float64, 1); conformity=:L2)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 4)

  lhs(u, v) = ∫(u * v) * dΩ
  rhs(v)    = ∫(v) * dΩ

  ch = ConstraintHandler(V, W, lhs, rhs)
  close!(ch)

  @test num_constrained_dofs(ch) == 2
  @test num_free_dofs(ch)        == 1

  x = zeros(ch.n_dofs)
  x[only(filter(i -> !is_constrained(ch, i), 1:ch.n_dofs))] = 2.0
  apply_constraints!(x, ch)
  uh = FEFunction(V, x)
  xf = CellField(p -> p[1], Ω)
  @test sum(∫(uh) * dΩ)       ≈ 1.0  atol=1e-12
  @test sum(∫(uh * xf) * dΩ) ≈ 0.5  atol=1e-12
end

let
  model = CartesianDiscreteModel((0,1), (3,))
  V = FESpace(model, ReferenceFE(lagrangian, Float64, 2); conformity=:L2)
  W = FESpace(model, ReferenceFE(lagrangian, Float64, 1); conformity=:L2)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 4)
  lhs(u, v) = ∫(u * v) * dΩ
  rhs(v)    = ∫(v) * dΩ
  ch = ConstraintHandler(V, W, lhs, rhs)
  close!(ch)
  @test num_constrained_dofs(ch) == 6
  @test num_free_dofs(ch)        == 3
end

# ===========================================================================
# Constructor 2 — cell, linear: cell-mean repair, f(x)=2x
# ===========================================================================

let
  model = CartesianDiscreteModel((0,1), (3,))
  V  = FESpace(model, ReferenceFE(lagrangian, Float64, 1); conformity=:L2)
  W0 = FESpace(model, ReferenceFE(lagrangian, Float64, 0); conformity=:L2)
  Ω  = Triangulation(model)
  dΩ = Measure(Ω, 4)

  f      = CellField(x -> 2 * x[1], Ω)
  lhs(v) = ∫(v) * dΩ
  rhs()  = ∫(f) * dΩ

  ch = ConstraintHandler(V, lhs, rhs)
  close!(ch)

  n_cells = num_cells(model)   # 3
  n_dofs  = ch.n_dofs          # 6

  @test num_constrained_dofs(ch) == n_cells
  @test num_free_dofs(ch)        == n_cells

  x = zeros(n_dofs)
  for (k, i) in enumerate(filter(i -> !is_constrained(ch, i), 1:n_dofs))
    x[i] = Float64(k)
  end
  apply_constraints!(x, ch)
  uh = FEFunction(V, x)

  for cell in 1:n_cells
    ind_vec = zeros(n_cells); ind_vec[cell] = 1.0
    ind = FEFunction(W0, ind_vec)
    @test sum(∫(uh * ind) * dΩ) ≈ sum(∫(f * ind) * dΩ)  atol=1e-12
  end
end

# ===========================================================================
# Constructor 3 — patch, bilinear: zero-mean per patch, 4×2 mesh
# ===========================================================================

let
  model = CartesianDiscreteModel((0,1,0,1), (4,2))
  V = FESpace(model, ReferenceFE(lagrangian, Float64, 1))

  topo       = Geometry.get_grid_topology(model)
  patch_data = Int32[1,2,5,6,  3,4,7,8]
  patch_ptrs = Int32[1,5,9]
  ptopo = Geometry.PatchTopology(topo, Table(patch_data, patch_ptrs))

  Ωp  = Geometry.PatchTriangulation(model, ptopo)
  dΩp = Measure(Ωp, 4)

  W_const   = ConstantFESpace(model)
  g         = CellField(x -> 0.0, Ωp)
  lhs(u, v) = ∫(u * v) * dΩp
  rhs(v)    = ∫(g * v) * dΩp

  ch = ConstraintHandler(ptopo, V, W_const, lhs, rhs)
  close!(ch)

  n_fdofs   = num_free_dofs(V)
  n_dofs    = n_fdofs + num_dirichlet_dofs(V)
  n_patches = Geometry.num_patches(ptopo)

  @test num_constrained_dofs(ch) == n_patches
  @test num_free_dofs(ch)        == n_fdofs - n_patches

  x = zeros(n_dofs)
  for i in 1:n_dofs
    !is_constrained(ch, i) && (x[i] = Float64(i))
  end
  apply_constraints!(x, ch)

  uh = FEFunction(V, x[1:n_fdofs])
  Ω  = Triangulation(model)
  dΩ = Measure(Ω, 4)
  W0 = FESpace(model, ReferenceFE(lagrangian, Float64, 0); conformity=:L2)
  patch_cells_tbl = Geometry.get_patch_cells(ptopo)
  for p in 1:n_patches
    cells_p       = collect(Int, patch_cells_tbl[p])
    indicator_vec = zeros(num_cells(model))
    indicator_vec[cells_p] .= 1.0
    indicator = FEFunction(W0, indicator_vec)
    @test sum(∫(uh * indicator) * dΩ) ≈ 0.0  atol=1e-10
  end
end

# ===========================================================================
# Constructor 4 — patch, linear: patch-mean = f(x,y)=x+y, 4×2 mesh
# ===========================================================================

let
  model = CartesianDiscreteModel((0,1,0,1), (4,2))
  V = FESpace(model, ReferenceFE(lagrangian, Float64, 1))

  topo       = Geometry.get_grid_topology(model)
  patch_data = Int32[1,2,5,6,  3,4,7,8]
  patch_ptrs = Int32[1,5,9]
  ptopo = Geometry.PatchTopology(topo, Table(patch_data, patch_ptrs))

  Ωp  = Geometry.PatchTriangulation(model, ptopo)
  dΩp = Measure(Ωp, 4)

  f      = CellField(x -> x[1] + x[2], Ωp)
  lhs(v) = ∫(v) * dΩp
  rhs()  = ∫(f) * dΩp

  ch = ConstraintHandler(ptopo, V, lhs, rhs)
  close!(ch)

  n_patches = Geometry.num_patches(ptopo)
  n_fdofs   = num_free_dofs(V)

  @test num_constrained_dofs(ch) == n_patches
  @test num_free_dofs(ch)        == n_fdofs - n_patches

  n_dofs = ch.n_dofs
  x = zeros(n_dofs)
  for i in 1:n_dofs
    !is_constrained(ch, i) && (x[i] = Float64(i))
  end
  apply_constraints!(x, ch)
  uh = FEFunction(V, x[1:n_fdofs])

  Ω  = Triangulation(model)
  dΩ = Measure(Ω, 4)
  W0 = FESpace(model, ReferenceFE(lagrangian, Float64, 0); conformity=:L2)
  g  = CellField(x -> x[1] + x[2], Ω)
  patch_cells_tbl = Geometry.get_patch_cells(ptopo)
  for p in 1:n_patches
    cells_p       = collect(Int, patch_cells_tbl[p])
    indicator_vec = zeros(num_cells(model))
    indicator_vec[cells_p] .= 1.0
    indicator = FEFunction(W0, indicator_vec)
    @test sum(∫(uh * indicator) * dΩ) ≈ sum(∫(g * indicator) * dΩ)  atol=1e-10
  end
end

end # module