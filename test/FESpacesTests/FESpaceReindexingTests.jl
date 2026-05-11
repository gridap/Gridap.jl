module FESpaceReindexingTests

using Test
using SparseArrays
using Gridap.Arrays
using Gridap.Geometry
using Gridap.ReferenceFEs
using Gridap.FESpaces
using Gridap.CellData

import Gridap.Arrays: compute_adjacency       # internal, accessed explicitly for adj kwarg test
import Gridap.FESpaces: compute_dof_permutation  # internal, accessed explicitly

model = CartesianDiscreteModel((0,1,0,1),(4,4))
reffe = ReferenceFE(lagrangian, Float64, 2)
V     = FESpace(model, reffe; dirichlet_tags="boundary")
nf    = num_free_dofs(V)
ndir  = num_dirichlet_dofs(V)

# reindex_free_dof_ids (algorithm dispatch) preserves DOF counts
for alg in (:rcm, :sloan, :coordinates)
  Vr = reindex_free_dof_ids(V, alg)
  @test num_free_dofs(Vr)      == nf
  @test num_dirichlet_dofs(Vr) == ndir
end

# :coordinates with custom by (sort by x-coordinate only)
Vx = reindex_free_dof_ids(V, :coordinates; by=x->x[1])
@test num_free_dofs(Vx) == nf

# reindex_free_dof_ids (direct permutation) preserves DOF counts
adj = compute_adjacency(get_cell_dof_ids(V), nf)
Vd  = reindex_free_dof_ids(V, compute_dof_permutation(adj, :rcm))
@test num_free_dofs(Vd)      == nf
@test num_dirichlet_dofs(Vd) == ndir

# Pre-built adj table can be passed explicitly
Vp = reindex_free_dof_ids(V, :rcm; adj=adj)
@test num_free_dofs(Vp) == nf

# RCM reduces bandwidth (end-to-end correctness check)
Ω  = Triangulation(model)
dΩ = Measure(Ω, 5)
a(u, v) = ∫(∇(v) ⋅ ∇(u))dΩ
A  = assemble_matrix(a, V, V)
Vr = reindex_free_dof_ids(V, :rcm)
Ar = assemble_matrix(a, Vr, Vr)
bw(M) = maximum(abs(i-j) for (i,j,_) in zip(findnz(M)...))
@test bw(Ar) < bw(A)

end # module
