module GatherAllocationRegressionTests

# Regression tests for the gather_free_values! / gather_dirichlet_values! fix.
#
# Before the fix, both functions allocated a throwaway vector that was
# immediately discarded, and this throwaway scaled with mesh size:
#
#   gather_free_values!      allocated a throwaway zero_dirichlet_values (~640–3200 B)
#   gather_dirichlet_values! allocated a throwaway zero_free_values    (~648–78400 B)
#
# Fix in src/FESpaces/UnconstrainedFESpaces.jl: added gather_free_values! and
# gather_dirichlet_values! overrides that call dedicated _gather_only_*_fill!
# helpers instead of delegating to gather_free_and_dirichlet_values!.
#
# After the fix, allocations are ~3500 bytes (constant, from array_cache setup)
# regardless of mesh size.

using Test
using Gridap
using Gridap.FESpaces
using Gridap.ReferenceFEs

function make_trial_space(n)
  model = CartesianDiscreteModel((0,1,0,1),(n,n))
  V = TestFESpace(model, ReferenceFE(lagrangian,Float64,1); dirichlet_tags="boundary")
  TrialFESpace(V, 0.0)
end

# ---------------------------------------------------------------------------
# Test 1: allocations no longer scale with free/dirichlet DOF count
# ---------------------------------------------------------------------------

@testset "gather allocation regression (20×20)" begin
  U  = make_trial_space(20)
  fv = zeros(num_free_dofs(U))
  dv = zeros(num_dirichlet_dofs(U))
  cv = get_cell_dof_values(zero(U))

  for _ in 1:3
    gather_free_values!(fv, U, cv)
    gather_dirichlet_values!(dv, U, cv)
  end

  alloc_free = @allocated gather_free_values!(fv, U, cv)
  alloc_dir  = @allocated gather_dirichlet_values!(dv, U, cv)

  println("gather_free_values!      : $alloc_free bytes")
  println("gather_dirichlet_values! : $alloc_dir bytes")

  # Before fix: ~4320 / ~6624 bytes (throwaway vectors included).
  # After fix:  ~3584 / ~3552 bytes (array_cache setup only, constant).
  # Threshold chosen above the unavoidable array_cache cost (~3500 B)
  # and below the pre-fix values.
  @test alloc_free < 5000
  @test alloc_dir  < 5000
end

# ---------------------------------------------------------------------------
# Test 2: allocations are constant across mesh sizes (the throwaway scaled)
# ---------------------------------------------------------------------------

@testset "gather allocation is mesh-size independent" begin
  println("\nMesh    | Free DOFs | Dir DOFs | gather_free (B) | gather_dir (B)")
  println("--------|-----------|----------|-----------------|---------------")

  results_free = Int[]
  results_dir  = Int[]

  for n in (10, 20, 50)
    U  = make_trial_space(n)
    fv = zeros(num_free_dofs(U))
    dv = zeros(num_dirichlet_dofs(U))
    cv = get_cell_dof_values(zero(U))

    for _ in 1:3
      gather_free_values!(fv, U, cv)
      gather_dirichlet_values!(dv, U, cv)
    end

    af = @allocated gather_free_values!(fv, U, cv)
    ad = @allocated gather_dirichlet_values!(dv, U, cv)
    push!(results_free, af)
    push!(results_dir,  ad)

    nf = num_free_dofs(U)
    nd = num_dirichlet_dofs(U)
    println("$(lpad(n,3))×$(lpad(n,3))   | $(lpad(nf,9)) | $(lpad(nd,8)) | $(lpad(af,15)) | $(lpad(ad,14))")

    @test af < 5000
    @test ad < 5000
  end

  # Before fix, gather_dirichlet_values! scaled with n_free (throwaway free vector).
  # After fix, both are constant. Check they don't grow more than 10% from 10×10 to 50×50.
  @test results_free[end] < 1.1 * results_free[1]
  @test results_dir[end]  < 1.1 * results_dir[1]
end

# ---------------------------------------------------------------------------
# Test 3: correctness – results match gather_free_and_dirichlet_values!
# ---------------------------------------------------------------------------

@testset "gather correctness" begin
  U  = make_trial_space(20)
  cv = get_cell_dof_values(zero(U))

  fv_ref = zero_free_values(U)
  dv_ref = zero_dirichlet_values(U)
  gather_free_and_dirichlet_values!(fv_ref, dv_ref, U, cv)

  fv = zeros(num_free_dofs(U))
  dv = zeros(num_dirichlet_dofs(U))
  gather_free_values!(fv, U, cv)
  gather_dirichlet_values!(dv, U, cv)

  @test fv == fv_ref
  @test dv == dv_ref

  # Non-trivial Dirichlet data
  model = CartesianDiscreteModel((0,1,0,1),(20,20))
  V2 = TestFESpace(model, ReferenceFE(lagrangian,Float64,1); dirichlet_tags="boundary")
  U2 = TrialFESpace(V2, x -> x[1] + x[2])
  uh = interpolate(x -> x[1] + x[2], U2)
  cv2 = get_cell_dof_values(uh)

  fv_ref2 = zero_free_values(U2)
  dv_ref2 = zero_dirichlet_values(U2)
  gather_free_and_dirichlet_values!(fv_ref2, dv_ref2, U2, cv2)

  fv2 = zeros(num_free_dofs(U2))
  dv2 = zeros(num_dirichlet_dofs(U2))
  gather_free_values!(fv2, U2, cv2)
  gather_dirichlet_values!(dv2, U2, cv2)

  @test fv2 ≈ fv_ref2
  @test dv2 ≈ dv_ref2
end

# ---------------------------------------------------------------------------
# Test 4: full interpolate! chain
# ---------------------------------------------------------------------------

@testset "interpolate! chain allocations" begin
  U  = make_trial_space(20)
  fv = zero_free_values(U)
  f  = x -> x[1]^2 + x[2]^2

  for _ in 1:3; interpolate!(f, fv, U); end

  alloc = @allocated interpolate!(f, fv, U)
  println("interpolate!(f, fv, U) : $alloc bytes  (was ~19648 before fix)")

  # Remaining allocations are from _cell_vals / CellField construction,
  # not the gather step. Threshold confirms the gather throwaway is gone.
  @test alloc < 19000
end

end # module
