module FormValuedPΛFESpacesTests

# Mesh level tests of the form valued (`ExteriorFormValue`) FEEC elements, i.e.
# of `compute_cell_bases_changes(::ReferenceFEName, ::TensorPushforward, ...)`
# in src/FESpaces/Pullbacks.jl.
#
# The reference level is covered by test/ReferenceFEsTests/FormValuedPΛRefFEsTests.jl;
# here everything is about what happens once the elements are glued over a mesh
# and transported to the physical cells.

using Test
using LinearAlgebra
using Random
using Gridap
using Gridap.FESpaces
using Gridap.Geometry
using Gridap.TensorValues
using Gridap.CellData
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Arrays

# The Hodge pairing ⟨ω,η⟩ = ★(ω ∧ ★η).
eucl(a,b) = inner(a,b,I)
L2norm(d,dΩ) = sqrt(abs(sum( ∫( Operation(eucl)(d,d) )dΩ )))

# The test forms below are not normalised, so every error is measured relative
# to the form it approximates.
relative_error(uh,ω,dΩ) = L2norm(uh - ω, dΩ) / L2norm(ω, dΩ)

# A degree `q` polynomial, `s` selecting one of a family of them.
_coef(x,q,s) = sum( (s*j+1) * sum((j+i+s)*x[i] for i in 1:length(x))^j for j in 0:q )

# A polynomial `k`-form in dimension `D` with coefficients of degree `q`. For
# `q` within reach of the space under test it lies in 𝓟_qΛᵏ ⊆ the space, so a
# conforming space with the right transport reproduces it exactly.
# `eltype(x)`, not `Float64`: the exact derivative is taken by differentiating
# this through the `Field` machinery, which evaluates it on dual numbers.
_form(::Val{k},::Val{D},q) where {k,D} =
  x -> ExteriorFormValue{k,D,eltype(x)}(ntuple(s -> _coef(x,q,s), binomial(D,k))...)

# The highest coefficient degree the space `𝓕_rΛᵏ` is guaranteed to contain:
# `r` for the full family, and only `r-1` for the trimmed one, which is complete
# one order below -- except at k=0, where 𝓟⁻_rΛ⁰ = 𝓟_rΛ⁰.
_full_degree(F,r,k) = (F == :P || k == 0) ? r : r-1

const TOL = 1.0e-12

sorted_tri = simplexify(CartesianDiscreteModel((0,1,0,1),(2,2)))
sorted_tet = simplexify(CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2)))

# Cells of five different areas and of both orientations, so that neither the
# `1/det(J)` factor nor its sign can hide behind a uniform mesh.
graded_tri = simplexify(CartesianDiscreteModel((0,1,0,1),(3,3); map = x -> VectorValue(x[1]^2,x[2])))

# ─────────────────────────────────────────────────────────────────────────────
# The global space is the one of the vector proxy: same number of free DOFs
# ─────────────────────────────────────────────────────────────────────────────
#
# This is the half of the correctness statement that the projection below cannot
# make: a space that is not glued at all, or glued too little, still contains the
# polynomial form and still projects it exactly -- it just has too many DOFs.

@testset "dof count parity with the vector proxies" begin
  for (model,p,D) in ((sorted_tri,TRI,2), (sorted_tet,TET,3)), r in 1:2
    proxies = (
      (:P , 0  , lagrangian     , r  ),
      (:P⁻, 1  , nedelec        , r-1),
      (:P , 1  , nedelec2       , r  ),
      (:P⁻, D-1, raviart_thomas , r-1),
      (:P , D-1, bdm            , r  ),
    )
    for (F,k,proxy_name,proxy_order) in proxies
      form  = FESpace(model, ReferenceFE(p,F,r,k,Float64; DG_calc=true))
      proxy = FESpace(model, ReferenceFE(p,proxy_name,Float64,proxy_order))
      @test num_free_dofs(form) == num_free_dofs(proxy)
    end
  end
end

# ─────────────────────────────────────────────────────────────────────────────
# A polynomial form of the space is reproduced exactly
# ─────────────────────────────────────────────────────────────────────────────
#
# This is the other half: it fails as soon as the space is glued *too much*, or
# the transport is wrong, because the global space then no longer contains the
# globally smooth polynomial form.

function projection_error(model,p,::Val{D},F,r,::Val{k}) where {D,k}
  q = _full_degree(F,r,k)
  V = FESpace(model, ReferenceFE(p,F,r,k,Float64; DG_calc=true))
  Ω = Triangulation(model)
  dΩ = Measure(Ω, 2*(r+q)+2)
  ω = CellField(_form(Val(k),Val(D),q), Ω)
  op = AffineFEOperator(
    (u,v) -> ∫( Operation(eucl)(u,v)  )dΩ,
     v    -> ∫( Operation(eucl)(ω,v)  )dΩ, V, V)
  relative_error(solve(op), ω, dΩ)
end

@testset "exact L2 projection on $lbl" for (lbl,model,p,D,rs) in (
    ("sorted TRI", sorted_tri, TRI, 2, 1:3),
    ("graded TRI", graded_tri, TRI, 2, 1:2),
    ("sorted TET", sorted_tet, TET, 3, 1:1),
  )
  for F in (:P,:P⁻), r in rs, k in 0:D
    @test projection_error(model,p,Val(D),F,r,Val(k)) < TOL
  end
end

# ─────────────────────────────────────────────────────────────────────────────
# The trimmed cubical family 𝓠⁻ᵣΛᵏ(□ᴰ)
# ─────────────────────────────────────────────────────────────────────────────

cart_quad   = CartesianDiscreteModel((0,1,0,1),(2,2))
cart_hex    = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))
graded_quad = CartesianDiscreteModel((0,1,0,1),(4,4); map = x -> VectorValue(x[1]^2, 0.5*x[2]))

@testset "Q⁻ dof count parity with the vector proxies" begin
  for (model,p,D) in ((cart_quad,QUAD,2), (cart_hex,HEX,3)), r in 1:2
    proxies = (
      (0  , lagrangian    , r  ),
      (1  , nedelec       , r-1),
      (D-1, raviart_thomas, r-1),
    )
    for (k,proxy_name,proxy_order) in proxies
      form  = FESpace(model, ReferenceFE(p,:Q⁻,r,k,Float64; DG_calc=true))
      proxy = FESpace(model, ReferenceFE(p,proxy_name,Float64,proxy_order))
      @test num_free_dofs(form) == num_free_dofs(proxy)
    end
  end
end

@testset "Q⁻ exact L2 projection on $lbl" for (lbl,model,p,D,rs) in (
    ("cartesian QUAD", cart_quad,   QUAD, 2, 1:3),
    ("graded QUAD",    graded_quad, QUAD, 2, 1:2),
    ("cartesian HEX",  cart_hex,    HEX,  3, 1:1),
  )
  for r in rs, k in 0:D
    @test projection_error(model,p,Val(D),:Q⁻,r,Val(k)) < TOL
  end
end

# ─────────────────────────────────────────────────────────────────────────────
# No `NormalSignMap`: the conforming trace is single valued as it stands
# ─────────────────────────────────────────────────────────────────────────────
# The AFW moments read `ι_f̂^* ω`, pulled back along the face's own parametrisation
# (`form_moment`), not `u⋅n̂` as the proxies do: no side is chosen, so no sign is
# owed, and there is no `apply_face_signflip` for `PΛConformity` either.

"component ids tangential (resp. not) to a facet of normal `eⱼ`, for each `j`"
_comp_ids(D,k,tangential) =
  [[s for (s,I) in enumerate(Gridap.TensorValues.sorted_combinations(D,k))
       if (j in I) != tangential] for j in 1:D]

_nrm_dir(n) = argmax(ntuple(i -> abs(n[i]), length(n)))

function _jump2(ωp, ωm, np, ids_per_dir)
  ids = ids_per_dir[_nrm_dir(np)]
  isempty(ids) && return zero(ωp.data[1])
  sum((ωp.data[s] - ωm.data[s])^2 for s in ids)
end

@testset "Q⁻Λᵏ needs no NormalSignMap on $lbl" for (lbl,model,p,D,krs) in (
    ("cartesian QUAD", cart_quad,   QUAD, 2, ((1,1),(1,2),(1,3))),
    ("graded QUAD",    graded_quad, QUAD, 2, ((1,1),(1,2))),
    ("cartesian HEX",  cart_hex,    HEX,  3, ((1,1),(1,2),(2,1),(2,2))),
  )
  for (k,r) in krs
    reffe = ReferenceFE(p,:Q⁻,r,k,Float64; DG_calc=true)
    V = FESpace(model, reffe)
    Random.seed!(1234)
    uh = FEFunction(V, rand(num_free_dofs(V)))

    Γ  = SkeletonTriangulation(model)
    dΓ = Measure(Γ, 2*(r+2))
    n  = get_normal_vector(Γ)
    J(u,ids) = sum(∫( Operation((a,b,c) -> _jump2(a,b,c,ids))(u.plus,u.minus,n.plus) )dΓ)

    @test J(uh, _comp_ids(D,k,true))  < 1.0e-20  # on `eⱼ` the trace keeps `dx^I`, `j ∉ I`
    @test J(uh, _comp_ids(D,k,false)) > 1.0e-6   # calibration: the dropped components jump

    # Vacuous where the sign map is the identity: it only touches facet owned DOFs.
    cell_reffe = [reffe for _ in 1:num_cells(model)]
    sgn = lazy_map(d -> diag(d), Gridap.FESpaces.get_sign_flip(model, cell_reffe))
    if any(c -> any(<(0), sgn[c]), 1:num_cells(model))
      basis = get_fe_basis(V)
      cvb = lazy_map(Broadcasting(*), get_cell_dof_values(uh), sgn)
      ub  = GenericCellField(lazy_map(linear_combination, cvb, get_data(basis)),
                             get_triangulation(basis), DomainStyle(basis))
      @test J(ub, _comp_ids(D,k,true)) > 1.0e-6  # calibration: a sign flip breaks it
    end
  end
end

@testset "exact interpolation" begin
  for F in (:P,:P⁻), r in 1:2, k in 0:2
    q = _full_degree(F,r,k)
    V = FESpace(sorted_tri, ReferenceFE(TRI,F,r,k,Float64; DG_calc=true))
    Ω = Triangulation(sorted_tri); dΩ = Measure(Ω, 2*(r+q)+2)
    f = _form(Val(k),Val(2),q)
    @test relative_error(interpolate(f,V), CellField(f,Ω), dΩ) < TOL
  end
end

# ─────────────────────────────────────────────────────────────────────────────
# Differentiating a form valued FE function
# ─────────────────────────────────────────────────────────────────────────────
# `∇` is refused by design: the `Broadcasting(∇)` reroute of
# src/ReferenceFEs/Pullbacks.jl hands `TensorPushforward` a reference gradient
# whose derivation index carries no variance, because `gradient_type` builds it
# out of a `Point`. `∂` is the spelling that works: its own reroute maps the
# value index block of the reference partials with `push_second_partials`, and
# `push_first_partials` transports their derivation index when the cell field
# is differentiated, so each index block is transported exactly once.
#
# The test forms lie in the space, so the FE derivative must reproduce the
# exact one to round-off. The errors are measured relative to the exact
# derivative, which is not normalised.

@testset "derivatives of a form valued FE function" begin
  cell_vals(cf,x) = (a = cf(x); c = array_cache(a); [getindex!(c,a,i) for i in 1:length(a)])
  _maxabs(vs) = maximum(maximum(abs.(get_array(v))) for cell in vs for v in cell)
  _maxdiff(a,b) = maximum(maximum(abs.(get_array(ai) .- get_array(bi)))
                          for (as,bs) in zip(a,b) for (ai,bi) in zip(as,bs))
  _err(a,b) = (s = _maxabs(b); s > 0 ? _maxdiff(a,b)/s : _maxdiff(a,b))

  for F in (:P,:P⁻), r in 1:2, k in 0:2
    q = _full_degree(F,r,k)
    V = FESpace(sorted_tri, ReferenceFE(TRI,F,r,k,Float64; DG_calc=true))
    Ω = Triangulation(sorted_tri); dΩ = Measure(Ω, 2*(r+q)+2)
    f = _form(Val(k),Val(2),q)
    uh = interpolate(f, V); fc = CellField(f, Ω)
    x = get_cell_points(dΩ)

    if k < 2
      @test_throws ArgumentError ∇(uh)(x)
    else
      # at k=D the reference gradient has no `pushforward` method at all, so
      # the variance check that refuses the lower degrees is never reached
      @test_throws MethodError ∇(uh)(x)
    end

    @test _err(cell_vals(∂(uh),x), cell_vals(∂(fc),x)) < TOL
    # 𝑑 of a 2-form vanishes in 2D, and the operator refuses that degree
    k < 2 && @test _err(cell_vals(𝑑(uh),x), cell_vals(𝑑(fc),x)) < TOL
  end
end

# `𝑑` of a cell basis takes the `natural_exterior_derivative` short circuit,
# which the FE function tests above do not reach in the same shape. On a 0-form
# space `∫ 𝑑u ∧ ⋆𝑑v` is the stiffness matrix `∫ ∇u⋅∇v` of the P1 Lagrangian
# space, up to the DOF numbering, which the spectrum is blind to.

@testset "FEEC stiffness matrix of a 0-form space" begin
  Ω = Triangulation(sorted_tri); dΩ = Measure(Ω, 4)
  V = FESpace(sorted_tri, ReferenceFE(TRI,:P,1,0,Float64; DG_calc=true))
  U = TrialFESpace(V)
  A = assemble_matrix((u,v) -> ∫( vol_coeff(𝑑(u) ∧ flat_hodge_star(𝑑(v))) )dΩ, U, V)

  Vl = FESpace(sorted_tri, ReferenceFE(lagrangian,Float64,1))
  Ul = TrialFESpace(Vl)
  Al = assemble_matrix((u,v) -> ∫( ∇(u)⋅∇(v) )dΩ, Ul, Vl)

  @test size(A) == size(Al)
  @test eigvals(Symmetric(Array(A))) ≈ eigvals(Symmetric(Array(Al)))
end

# ─────────────────────────────────────────────────────────────────────────────
# The mass matrix of the Hodge pairing
# ─────────────────────────────────────────────────────────────────────────────

@testset "Hodge mass matrix is symmetric positive definite" begin
  Ω = Triangulation(sorted_tri); dΩ = Measure(Ω,6)
  for F in (:P,:P⁻), r in 1:2, k in 0:2
    V = FESpace(sorted_tri, ReferenceFE(TRI,F,r,k,Float64; DG_calc=true))
    u, v = get_trial_fe_basis(V), get_fe_basis(V)
    M = Matrix(assemble_matrix(∫( Operation(eucl)(u,v) )dΩ, V, V))
    @test M ≈ M'
    @test isposdef(Symmetric(M))
  end
end

# ─────────────────────────────────────────────────────────────────────────────
# The `D`-form element -- the "DG" end of the table
# ─────────────────────────────────────────────────────────────────────────────

@testset "D-form element" begin
  V = FESpace(sorted_tri, ReferenceFE(TRI,:P⁻,1,2,Float64; DG_calc=true))
  # 𝓟⁻₁Λ² is spanned by the piecewise constant densities: one DOF per cell, and
  # none of them glued -- the space is discontinuous.
  @test num_free_dofs(V) == num_cells(sorted_tri)

  # The single DOF of a cell is a moment ∫_K̂ (ι*ω) ∧ μ against one fixed reference
  # 0-form μ, so by duality ∫_K φ_K is the *same* number on every cell however
  # large the cell is -- only its sign changes, with the cell orientation, as the
  # integral of a form must. That only holds because the `1/det(J)` transport is
  # applied; with the identity transport ∫_K φ_K would be proportional to |K|.
  # (The number itself depends on how μ is normalised, so it is not asserted.)
  Ω = Triangulation(graded_tri); dΩ = Measure(Ω,4)
  W = FESpace(graded_tri, ReferenceFE(TRI,:P⁻,1,2,Float64; DG_calc=true))
  areas = get_array(∫(1)dΩ)
  @test maximum(areas)/minimum(areas) > 4   # the mesh really is graded

  uh = FEFunction(W, ones(num_free_dofs(W)))
  cell_integrals = get_array(∫( Operation(ω -> first(ω.data))(uh) )dΩ)
  @test all(c -> abs(abs(c) - abs(first(cell_integrals))) < TOL, cell_integrals)
end

# ─────────────────────────────────────────────────────────────────────────────
# The nodal names: the two ends of the table, on simplices and on n-cubes
# ─────────────────────────────────────────────────────────────────────────────
#
# `k = 0` and `k = D` are the only degrees an n-cube reaches so far (the
# `0 < k < D` ones need the `:Q⁻` form valued bases), and the only ones a nodal
# element is defined for.

quad_model = CartesianDiscreteModel((0,1,0,1),(2,2))
hex_model  = CartesianDiscreteModel((0,1,0,1,0,1),(2,2,2))

@testset "nodal form valued elements" begin
  for (model,p,D,F,r,k) in (
      (sorted_tri,TRI,2,:P ,2,0), (sorted_tri,TRI,2,:P ,2,2),
      (sorted_tri,TRI,2,:P ,3,2), (sorted_tet,TET,3,:P ,2,3),
      (quad_model,QUAD,2,:Q⁻,2,0), (quad_model,QUAD,2,:Q⁻,2,2),
      (quad_model,QUAD,2,:S ,2,0), (quad_model,QUAD,2,:S ,2,2),
      (hex_model ,HEX ,3,:Q⁻,1,0), (hex_model ,HEX ,3,:Q⁻,2,3),
      (hex_model ,HEX ,3,:S ,2,0), (hex_model ,HEX ,3,:S ,2,3),
    )
    reffe = ReferenceFE(p,F,r,k,Float64; DG_calc=true, nodal=true)
    V = FESpace(model, reffe)

    # 𝓠⁻ᵣΛᴰ drops to order r-1; the others keep r. 𝓢ᵣΛᴰ ≡ 𝓟ᵣ.
    q = (k == D && F == :Q⁻) ? r-1 : r
    Ω = Triangulation(model); dΩ = Measure(Ω, 2*(r+q)+2)
    f = _form(Val(k),Val(D),q)
    @test relative_error(interpolate(f,V), CellField(f,Ω), dΩ) < TOL

    if k == D
      # A D-form owns every DOF on its cell: the space is fully discontinuous.
      @test num_free_dofs(V) == num_cells(model)*num_dofs(reffe)
    end
  end

  # The modal names on an n-cube go through the `:Q⁻` moment table.
  @test num_dofs(ReferenceFE(QUAD,:Q⁻,2,0,Float64; DG_calc=true)) == 9
end

# ─────────────────────────────────────────────────────────────────────────────
# Conformity symbols and Dirichlet boundary conditions
# ─────────────────────────────────────────────────────────────────────────────

@testset "conformity symbols" begin
  reffe = ReferenceFE(TRI,:P,2,1,Float64; DG_calc=true)
  @test Conformity(reffe) isa PΛConformity
  for sym in (:HΛ, :HLambda)
    @test num_free_dofs(FESpace(sorted_tri,reffe; conformity=sym)) == 72
  end
  # `:L2` drops the inter-cell gluing but NOT the transport: a form is still a
  # form. The pushforward is therefore chosen from the value type, not from the
  # (name, conformity) pair, which would hand these elements an
  # `IdentityPiolaMap` and leave their values in the reference domain.
  q = _full_degree(:P,2,1)
  Ω = Triangulation(sorted_tri); dΩ = Measure(Ω, 2*(2+q)+2)
  f = _form(Val(1),Val(2),q)
  for conf in (:L2, L2Conformity())
    V = FESpace(sorted_tri, reffe; conformity=conf)
    # fully discontinuous: one copy of the element per cell, nothing shared
    @test num_free_dofs(V) == num_cells(sorted_tri)*num_dofs(reffe)
    # and still transported, so a polynomial form is reproduced exactly
    @test relative_error(interpolate(f,V), CellField(f,Ω), dΩ) < TOL
  end
end

@testset "DOF scaling" begin
  # `get_dofscale_setter_function(reffe,::TensorPushforward)` scales a `k`-form DOF by
  # `hᵏ`. That is a change of basis applied to the DOFs and to the shape functions alike,
  # so it must leave the space -- and hence the interpolant -- untouched.
  for F in (:P,:P⁻), r in 1:2, k in 0:2
    q = _full_degree(F,r,k)
    reffe = ReferenceFE(TRI,F,r,k,Float64; DG_calc=true)
    V = FESpace(sorted_tri, reffe; scale_dof=true, global_meshsize=0.5)
    Ω = Triangulation(sorted_tri); dΩ = Measure(Ω, 2*(r+q)+2)
    f = _form(Val(k),Val(2),q)
    @test relative_error(interpolate(f,V), CellField(f,Ω), dΩ) < TOL
  end

  # Without a `global_meshsize` the local one is the `d`-root of a `d`-volume, which is
  # not available at a vertex -- so only the elements with no vertex owned DOF can do it.
  @test num_free_dofs(FESpace(sorted_tri, ReferenceFE(TRI,:P,2,1,Float64; DG_calc=true);
                              scale_dof=true)) == 72
  @test_throws ErrorException FESpace(sorted_tri, ReferenceFE(TRI,:P,2,0,Float64; DG_calc=true);
                                      scale_dof=true)
end

@testset "Dirichlet boundary conditions" begin
  for (k,ntotal) in ((0,25),(1,72))
    reffe = ReferenceFE(TRI,:P,2,k,Float64; DG_calc=true)
    V = FESpace(sorted_tri, reffe; dirichlet_tags="boundary")
    @test num_free_dofs(V) + num_dirichlet_dofs(V) == ntotal
    @test num_dirichlet_dofs(V) > 0
  end
end

# ─────────────────────────────────────────────────────────────────────────────
# Meshes whose cells do not list their vertices in increasing global id order
# ─────────────────────────────────────────────────────────────────────────────

@testset "unsorted cells" begin
  grid = get_grid(sorted_tri)
  relabel(f) = UnstructuredDiscreteModel(UnstructuredGrid(
    get_node_coordinates(grid),
    Table([f(i,collect(c)) for (i,c) in enumerate(get_cell_node_ids(grid))]),
    get_reffes(grid), get_cell_type(grid)))

  # Reversing *every* cell leaves the mesh consistently oriented -- each shared
  # face is still seen in the same vertex order by both of its cells -- it only
  # stops being sorted. Reversing every other cell is what actually breaks the
  # agreement, and is the mesh the gluing has to be tested against.
  reversed = relabel((i,c) -> reverse(c))
  mixed    = relabel((i,c) -> isodd(i) ? reverse(c) : c)

  "number of faces whose two cells list its vertices in a different order"
  function face_order_clashes(model)
    D = num_cell_dims(model); topo = get_grid_topology(model)
    cell_verts  = get_faces(topo,D,0)
    cell_facets = get_faces(topo,D,D-1)
    lface_lverts = get_faces(TRI,D-1,0)
    seen = Dict{Int,Vector{Int}}(); clashes = 0
    for cell in 1:num_cells(model), (lface,face) in enumerate(cell_facets[cell])
      order = [cell_verts[cell][j] for j in lface_lverts[lface]]
      haskey(seen,face) ? (seen[face] != order && (clashes += 1)) : (seen[face] = order)
    end
    clashes
  end
  @test face_order_clashes(sorted_tri) == 0
  @test face_order_clashes(reversed)   == 0
  @test face_order_clashes(mixed)       > 0

  @test Gridap.Geometry.has_sorted_cells(sorted_tri)
  @test Gridap.Geometry.has_sorted_cells(sorted_tet)
  @test Gridap.Geometry.has_sorted_cells(graded_tri)
  @test !Gridap.Geometry.has_sorted_cells(reversed)
  @test !Gridap.Geometry.has_sorted_cells(mixed)

  # The nodal elements declare a DOF permutation per vertex permutation, so they
  # are glued on any mesh, exactly as their scalar counterparts are.
  nodal_k0 = ReferenceFE(TRI,:P,3,0,Float64; DG_calc=true, nodal=true)
  @test Gridap.ReferenceFEs.has_complete_dof_permutations(nodal_k0)
  for model in (reversed, mixed)
    V = FESpace(model, nodal_k0)
    @test num_free_dofs(V) == 49
    Ω = Triangulation(model); dΩ = Measure(Ω,10)
    f = _form(Val(0),Val(2),3)
    @test relative_error(interpolate(f,V), CellField(f,Ω), dΩ) < TOL
  end

  # The moment based ones declare the identity alone. That includes the modal
  # `k = 0` element -- it is *not* the nodal one, and from r = 3 on it carries
  # two orientation sensitive moments per edge -- and the `D`-forms, which glue
  # nothing but are still read by the gluing code face by face.
  for reffe in (ReferenceFE(TRI,:P,3,0,Float64; DG_calc=true),
                ReferenceFE(TRI,:P,2,1,Float64; DG_calc=true),
                ReferenceFE(TRI,:P,2,2,Float64; DG_calc=true))
    @test !Gridap.ReferenceFEs.has_complete_dof_permutations(reffe)
    @test_throws ErrorException FESpace(mixed, reffe)
  end
end

@testset "nodal D-forms on an unsorted n-cube mesh" begin
  # An n-cube mesh cannot be `simplexify`d into the shape the moment based elements need,
  # so the nodal `D`-form -- which glues nothing at all -- must not depend on the vertex
  # order. It declares one empty permutation per vertex permutation on every glued face
  # exactly so that it does not, see `_top_form_face_own_nodes_permutations`.
  cart = UnstructuredDiscreteModel(CartesianDiscreteModel((0,1,0,1),(2,2)))
  grid = get_grid(cart)
  permuted = UnstructuredDiscreteModel(UnstructuredGrid(
    get_node_coordinates(grid),
    Table([isodd(i) ? collect(c)[[4,3,2,1]] : collect(c)
           for (i,c) in enumerate(get_cell_node_ids(grid))]),
    get_reffes(grid), collect(get_cell_type(grid))))
  @test !Gridap.Geometry.has_sorted_cells(permuted)

  reffe = ReferenceFE(QUAD,:Q⁻,2,2,Float64; DG_calc=true, nodal=true)
  @test Gridap.ReferenceFEs.has_complete_dof_permutations(reffe)

  V = FESpace(permuted, reffe)
  @test num_free_dofs(V) == num_cells(permuted)*num_dofs(reffe)
  Ω = Triangulation(permuted); dΩ = Measure(Ω,8)
  f = _form(Val(2),Val(2),1)
  @test relative_error(interpolate(f,V), CellField(f,Ω), dΩ) < TOL

  # Moment based: identity permutation only, so the `has_sorted_cells` guard refuses it.
  for k in 0:2
    mom = ReferenceFE(QUAD,:Q⁻,2,k,Float64; DG_calc=true)
    @test !Gridap.ReferenceFEs.has_complete_dof_permutations(mom)
    @test_throws ErrorException FESpace(permuted, mom)
  end
end

end # module
