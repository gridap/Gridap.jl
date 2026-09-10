module CartProdRefFEsTests

using LinearAlgebra
using Test
using FillArrays

using Gridap
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.FESpaces
using Gridap.CellData

############################################################################################
# Helpers
############################################################################################

unit_square(n) = simplexify(CartesianDiscreteModel((0, 1, 0, 1), (n, n)))

pts = [Point(0.3, 0.2), Point(0.1, 0.5), Point(0.25, 0.25)]

# the atoms: a scalar one whose DoFs differentiate twice, a scalar one that does
# not, and two vector-valued ones
morley = MorleyRefFE(Float64, TRI)
argyris = ArgyrisRefFE(Float64, TRI)
lag2 = ReferenceFE(TRI, lagrangian, Float64, 2)
lagv = ReferenceFE(TRI, lagrangian, VectorValue{2,Float64}, 1)
rt0 = ReferenceFE(TRI, raviart_thomas, Float64, 0)
rt1 = ReferenceFE(TRI, raviart_thomas, Float64, 1)
mtw2 = MardalTaiWintherRefFE(Float64, TRI)

# K = 3 on a 2D polytope everywhere it matters: an implementation that prepends
# the copy index instead of appending it transposes every gradient, and that is
# invisible when K == D.

############################################################################################
# The reference element
############################################################################################

function test_cp_reffe(atom, K)
  r = CartProdRefFE(atom, Val(K))
  n = num_dofs(atom)

  @test r isa GenericRefFE
  @test num_dofs(r) == K * n
  @test length(get_shapefuns(r)) == K * n
  @test length(get_prebasis(r)) == K * length(get_prebasis(atom))
  @test get_polytope(r) == get_polytope(atom)
  @test Conformity(r) == Conformity(atom)

  # duality: the whole point of the construction
  M = evaluate(get_dof_basis(r), get_shapefuns(r))
  @test M ≈ Matrix(I, K*n, K*n)

  # ownership is the base's, expanded copy-fastest
  own_a = get_face_own_dofs(atom)
  own_r = get_face_own_dofs(r)
  @test length(own_r) == length(own_a)
  for (oa, orr) in zip(own_a, own_r)
    @test orr == Int[K*(d-1)+c for d in oa for c in 1:K]
  end
  @test sort(vcat(own_r...)) == collect(1:K*n)
end

test_cp_reffe(morley, 2)
test_cp_reffe(morley, 3)
test_cp_reffe(argyris, 3)
test_cp_reffe(rt0, 3)
test_cp_reffe(rt1, 2)
test_cp_reffe(mtw2, 2)
test_cp_reffe(lag2, 3)
test_cp_reffe(lagv, 2)

@test_throws ErrorException CartProdRefFE(morley, Val(0))

############################################################################################
# Stacking a nodal DoF basis
#
# The result is another `LagrangianDofBasis` on the same nodes. Its layout is
# built explicitly rather than by `LagrangianDofBasis(W, nodes)`, which numbers
# with the node fastest where everything here runs the copy fastest, so what is
# checked is the layout: copy `c` of DoF `a` is at `K*(a-1)+c` and reads
# component `j + d*(c-1)`, with `j` the atom's component and `d` its count.
############################################################################################

function test_cp_nodal_dofs(atom, K)
  adb = get_dof_basis(atom)
  db = ReferenceFEs._cp_stack_dofs(adb, Val(K))
  d = num_indep_components(testitem(evaluate(get_shapefuns(atom), pts)))

  @test db isa LagrangianDofBasis
  @test get_nodes(db) == get_nodes(adb)
  for a in 1:num_dofs(atom), c in 1:K
    @test db.dof_to_node[K*(a-1)+c] == adb.dof_to_node[a]
    @test db.dof_to_comp[K*(a-1)+c] == adb.dof_to_comp[a] + d*(c-1)
  end
end

test_cp_nodal_dofs(lag2, 3)
test_cp_nodal_dofs(lagv, 2)

# the stacked space is Gridap's own vector-valued Lagrangian space: the two
# number their DoFs differently, but they interpolate identically
let K = 3, r = 2, model = unit_square(3)
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 2*r + 2)
  Vs = FESpace(model, CartProdRefFE(ReferenceFE(TRI, lagrangian, Float64, r), Val(K)))
  Vg = FESpace(model, ReferenceFE(lagrangian, VectorValue{K,Float64}, r))
  @test num_free_dofs(Vs) == num_free_dofs(Vg)

  u(x) = VectorValue(1.0 + x[1]^2, x[1]*x[2] - 2.0, 3.0 - x[2]^2)
  e = interpolate(u, Vs) - interpolate(u, Vg)
  @test sqrt(sum( ∫( e ⋅ e )dΩ )) < 1e-13
end

############################################################################################
# The stacked basis is the atom's, scattered copy-fastest
#
# This is the ordering convention, checked directly rather than through a
# consequence of it.
############################################################################################

function test_cp_basis_values(atom, K)
  b = get_shapefuns(atom)
  sb = ReferenceFEs.CartProdBasis{K}(b)
  E = representatives_of_componentbasis_dual(VectorValue{K,Float64})
  va = evaluate(b, pts)
  vs = evaluate(sb, pts)
  @test size(vs) == (length(pts), K * size(va, 2))
  for i in axes(va, 1), j in axes(va, 2), c in 1:K
    @test vs[i, K*(j-1)+c] ≈ outer(va[i, j], E[c])
  end
end

test_cp_basis_values(morley, 3)
test_cp_basis_values(rt0, 3)

# ... and the gradient obeys the *same* rule, which is the point of this variant:
# appending commutes with Gridap's prepending of derivative indices. For a scalar
# atom, ∇φ is a VectorValue{D} and the stacked gradient a TensorValue{D,K} with
# column c holding ∇φ -- not row c, which is what prepending would give.
function test_cp_basis_gradients(atom, K)
  b = get_shapefuns(atom)
  E = representatives_of_componentbasis_dual(VectorValue{K,Float64})
  ga = evaluate(Broadcasting(∇)(b), pts)
  gs = evaluate(Broadcasting(∇)(ReferenceFEs.CartProdBasis{K}(b)), pts)
  for i in axes(ga, 1), j in axes(ga, 2), c in 1:K
    g = gs[i, K*(j-1)+c]
    @test g ≈ outer(ga[i, j], E[c])
  end
end

test_cp_basis_gradients(morley, 3)
test_cp_basis_gradients(argyris, 3)

# second derivatives too, which is what the Argyris DoFs need
let b = get_shapefuns(argyris), K = 3
  E = representatives_of_componentbasis_dual(VectorValue{K,Float64})
  Ha = evaluate(Broadcasting(∇∇)(b), pts)
  Hs = evaluate(Broadcasting(∇∇)(ReferenceFEs.CartProdBasis{K}(b)), pts)
  for i in axes(Ha, 1), j in axes(Ha, 2), c in 1:K
    @test Hs[i, K*(j-1)+c] ≈ outer(Ha[i, j], E[c])
  end
end

# a vector atom's gradient is the case that separates the two variants: here it
# is still a plain append, (D,d) -> (D,d,K), where the row layout needs a middle
# insertion into (D,K,d).
let b = get_shapefuns(rt1), K = 3
  E = representatives_of_componentbasis_dual(VectorValue{K,Float64})
  ga = evaluate(Broadcasting(∇)(b), pts)
  gs = evaluate(Broadcasting(∇)(ReferenceFEs.CartProdBasis{K}(b)), pts)
  for i in axes(ga, 1), j in axes(ga, 2), c in 1:K
    g = gs[i, K*(j-1)+c]
    @test g ≈ outer(ga[i, j], E[c])
    for k in 1:2, jj in 1:2, cc in 1:K
      @test g[k, jj, cc] ≈ (cc == c ? ga[i, j][k, jj] : 0.0)
    end
  end
end

############################################################################################
# `tr` of the stacked gradient is the per-copy divergence
#
# The property the column layout buys. Gridap's `tr` on a third-order tensor
# traces the first two indices, so on (k,j,c) it contracts the derivative index
# against the atom's own -- one divergence per copy. Under the row layout the
# same trace hits the copy index instead, and for K ≠ D does not even typecheck.
############################################################################################

let b = get_shapefuns(rt1), K = 3
  E = representatives_of_componentbasis_dual(VectorValue{K,Float64})
  ga = evaluate(Broadcasting(∇)(b), pts)
  gs = evaluate(Broadcasting(∇)(ReferenceFEs.CartProdBasis{K}(b)), pts)
  for i in axes(ga, 1), j in axes(ga, 2), c in 1:K
    @test tr(gs[i, K*(j-1)+c]) ≈ tr(ga[i, j]) * E[c]
  end
end

############################################################################################
# The stacked DoFs are ordinary Gridap DoF bases
#
# Nothing is wrapped: slicing the field is absorbed into the moment weights.
############################################################################################

@test ReferenceFEs._cp_stack_dofs(get_dof_basis(mtw2), Val(2)) isa MomentBasedDofBasis
@test ReferenceFEs._cp_stack_dofs(get_dof_basis(morley), Val(2)) isa ReferenceFEs.ConcatenatedDofVector
@test ReferenceFEs._cp_stack_dofs(get_dof_basis(rt0), Val(2)) isa ReferenceFEs.LinearCombinationDofVector

# σ_{i,c}(u) = σ_i(π_c u): applying the stacked DoFs to a field that is copy c of
# an atom field reproduces the atom's DoFs in the rows of copy c, and zero
# elsewhere.
function test_cp_dofs_slice(atom, K)
  db = get_dof_basis(atom)
  sdb = ReferenceFEs._cp_stack_dofs(db, Val(K))
  b = get_shapefuns(atom)
  va = evaluate(db, b)                       # (n, n), the identity
  for c in 1:K
    inj = ReferenceFEs.CartProdBasis{K}(b)                # all copies at once
    vs = evaluate(sdb, inj)
    for i in axes(va, 1), j in axes(va, 2), cc in 1:K
      @test vs[K*(i-1)+c, K*(j-1)+cc] ≈ (cc == c ? va[i, j] : 0.0)
    end
    break                                    # the loop above already covers all c
  end
end

test_cp_dofs_slice(morley, 3)
test_cp_dofs_slice(rt0, 2)

############################################################################################
# The change of basis is kron(P, I_K)
############################################################################################

let K = 3, model = unit_square(3)
  atom = morley
  r = CartProdRefFE(atom, Val(K))
  cell_map = get_cell_map(get_grid(model))
  cell_Jt = lazy_map(Broadcasting(∇), cell_map)

  base_ch = compute_cell_bases_changes(
    get_name(atom), Pushforward(get_name(atom), Conformity(atom)),
    model, Fill(atom, num_cells(model)), cell_Jt)
  cp_ch = compute_cell_bases_changes(
    get_name(r), Pushforward(get_name(r), Conformity(r)),
    model, Fill(r, num_cells(model)), cell_Jt)
  @test !isnothing(base_ch) && !isnothing(cp_ch)
  P, Pit = base_ch
  Q, Qit = cp_ch
  for cell in 1:num_cells(model)
    @test Matrix(Q[cell]) ≈ kron(Matrix(P[cell]), Matrix(1.0I, K, K))
    @test Matrix(Qit[cell]) ≈ kron(Matrix(Pit[cell]), Matrix(1.0I, K, K))
  end
end

############################################################################################
# FE spaces: vector-valued Morley
#
# The scalar-atom application. A global quadratic in each component must be
# interpolated exactly -- this is the path that goes through `interpolate`, i.e.
# through the DoF basis acting on an arbitrary `CellField`.
############################################################################################

let K = 3, model = unit_square(4)
  V = FESpace(model, CartProdRefFE(morley, Val(K)))
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 6)

  topo = get_grid_topology(model)
  @test num_free_dofs(V) == K * (num_faces(topo, 0) + num_faces(topo, 1))

  u(x) = VectorValue(1.0 + 2*x[1] - x[2] + 3*x[1]^2,
                     3.0 - x[1] + 0.5*x[2]^2,
                     x[1]*x[2] - 2.0)
  uh = interpolate(u, V)
  e = uh - u
  @test sqrt(sum( ∫( e ⋅ e )dΩ )) < 1e-12
end

############################################################################################
# FE spaces: column-wise H(div) tensors
#
# The vector-atom application: K stacked Raviart--Thomas columns, the stress
# space of elasticity with weakly imposed symmetry. Each column is
# H(div)-conforming, so the traction n⋅σ is continuous while the full tensor
# jumps.
############################################################################################

let K = 2, model = unit_square(4)
  V = FESpace(model, CartProdRefFE(rt0, Val(K)))
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 6)

  topo = get_grid_topology(model)
  @test num_free_dofs(V) == K * num_faces(topo, 1)

  # Each *column* must be an RT0 field, i.e. of the form a + b*x with a single b.
  # `TensorValue` is column major, so the entries below are (m11, m21, m12, m22):
  # column 1 = (1 + x, 2 + y), column 2 = (-1 + 2x, 0.5 + 2y).
  σ(x) = TensorValue{2,2,Float64}(1.0 + x[1], 2.0 + x[2], -1.0 + 2*x[1], 0.5 + 2*x[2])
  σh = interpolate(σ, V)
  e = σh - σ
  @test sqrt(sum( ∫( e ⊙ e )dΩ )) < 1e-12

  # and Gridap's own `divergence` is the vector of the columns' divergences:
  # (1 + 1, 2 + 2). This is the call that is meaningless under the row layout.
  d = divergence(σh) - CellField(x -> VectorValue(2.0, 4.0), Ω)
  @test sqrt(sum( ∫( d ⋅ d )dΩ )) < 1e-10

  # a field outside the space, to see the conformity
  w(x) = TensorValue{2,2,Float64}(sin(2*x[1]), cos(3*x[2]), sin(x[1]+x[2]), x[1]^2)
  wh = interpolate(w, V)

  Λ = SkeletonTriangulation(model)
  dΛ = Measure(Λ, 6)
  n = get_normal_vector(Λ).⁺

  @test sqrt(sum( ∫( (n ⋅ jump(wh)) ⋅ (n ⋅ jump(wh)) )dΛ )) < 1e-12
  @test sqrt(sum( ∫( jump(wh) ⊙ jump(wh) )dΛ )) > 1e-3
end

end # module
