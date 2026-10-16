module FormValuedPΛRefFEsTests
# Form valued (ExteriorFormValue) FEEC reference elements PᵣΛᵏ and P⁻ᵣΛᵏ:
# the TensorPushforward and its reduction to the classical Piola maps, the AFW
# moment DOFs, and the reference level validation of the whole element.
#
# These elements have no name of their own: they are the Nedelec / Raviart-
# Thomas / BDM / ModalScalar elements, built with an `ExteriorFormValue` value
# type instead of a vector proxy one.

using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.ReferenceFEs
using Gridap.Fields
using Gridap.Arrays
using LinearAlgebra
using FillArrays: Fill
using Test

simplex(D) = D == 2 ? TRI : TET

# The form valued FᵣΛᵏ(p), through the FEEC factory.
pλ(p,r,k,F) = ReferenceFE(p, F, r, k, Float64; DG_calc=true)

# ─── TensorPushforward: reduction to the classical Piola maps ────────────────
# A single map must reproduce every classical Piola map, selected by the
# variance of the value alone.

@testset "TensorPushforward vs the classical Piola maps" begin
  tpf  = TensorPushforward()
  itpf = inverse_map(tpf)
  # `Jt` is what `∇` of a Lagrangian map gives: no variance. `Jtv` is what
  # `partial_derivatives` gives, and it is the only one `TensorPushforward`
  # accepts -- the classical Piola maps below keep taking the un-varianced one.
  Jt  = TensorValue{3,3}(2.0,0.3,0.1, 0.2,1.5,0.4, 0.05,0.2,3.0)
  Jtv = TensorValue{Cov(3),Con(3)}(Jt.data...)

  # an un-varianced jacobian is refused rather than silently re-tagged
  for k in 1:3
    ω = ExteriorFormValue{k,Cov(3)}(ntuple(i->Float64(i), binomial(3,k))...)
    @test_throws ArgumentError evaluate!(nothing, tpf, ω, Jt)
  end

  # k=0: a 0-form has no index to transport, so it is carried unchanged -- but
  # it is a form, not a `Real`, and the pullback re-tags its dimension, which
  # is what lets it be traced onto a face like any other degree
  z0 = ExteriorFormValue{0,Cov(3)}(2.5)
  @test evaluate!(nothing, tpf, z0, Jtv) == z0
  @test pullback(z0, TensorValue{Cov(2),Con(3)}(1.0,0.0, 0.0,1.0, 0.5,0.5)) ==
        ExteriorFormValue{0,Cov(2)}(2.5)
  @test pullback(z0, TensorValue{Cov(0),Con(3),Float64,0}()) == ExteriorFormValue{0,Cov(0)}(2.5)
  # ... and a bare `Real` is *not* a 0-form: it carries neither a dimension nor
  # a variance, so it cannot be traced onto a face nor take a Hodge star. It is
  # refused rather than passed through.
  @test_throws ErrorException evaluate!(nothing, tpf, 2.5, Jtv)

  # k=1: exactly the covariant Piola map
  ω1 = ExteriorFormValue{1,Cov(3)}(1.0,-2.0,0.5)
  cov = evaluate!(nothing, CoVariantPiolaMap(), VectorValue(1.0,-2.0,0.5), Jt)
  @test collect(evaluate!(nothing, tpf, ω1, Jtv).data) ≈ collect(cov.data)

  # k=D: multiplication by 1/det(J)
  @test first(evaluate!(nothing, tpf, ExteriorFormValue{3,Cov(3)}(1.0), Jtv).data) ≈ 1/det(Jt)

  # k=D-1: the contravariant Piola map read through the Hodge proxy, but with
  # a SIGNED determinant instead of meas(Jt)
  ω2 = ExteriorFormValue{2,Cov(3)}(1.0,-2.0,0.5)
  h  = flat_hodge_star(ω2)
  contra_signed = VectorValue(h.data...) ⋅ ((1/det(Jt)) * Jt)
  @test collect(flat_hodge_star(evaluate!(nothing, tpf, ω2, Jtv)).data) ≈
        collect(contra_signed.data)

  # the inverse map is a left inverse, for every form degree
  for k in 0:3
    ω = ExteriorFormValue{k,Cov(3)}(ntuple(i->Float64(i), binomial(3,k))...)
    rt = evaluate!(nothing, itpf, evaluate!(nothing, tpf, ω, Jtv), Jtv)
    @test collect(rt.data) ≈ collect(ω.data)
  end
end

# ─── Half varianced values are refused, not guessed ──────────────────────────
# `pushforward`/`pullback` transport EVERY index of their argument. A value
# whose variance is only partly declared cannot be transported without guessing
# how the remaining indices transform, so it is rejected.
#
# This is what `∇` of a form-valued field produces: `gradient_type` builds the
# derivation index out of a `Point`, which carries no variance.
# `partial_derivatives` is the operator that declares it `Cov`.

@testset "TensorPushforward refuses a partly varianced value" begin
  tpf = TensorPushforward()
  Jtv = TensorValue{Cov(3),Con(3)}(2.0,0.3,0.1, 0.2,1.5,0.4, 0.05,0.2,3.0)

  # `∇` of a 1-form valued field: index 1 is the derivation index (plain Int),
  # index 2 is the covariant value index
  g = TensorValue{3,Cov(3),Float64,9}(ntuple(i->Float64(i),9)...)
  @test_throws ArgumentError evaluate!(nothing, tpf, g, Jtv)
  @test_throws ArgumentError evaluate!(nothing, inverse_map(tpf), g, Jtv)

  # `∇` of a 0-form: no varianced index at all
  v = VectorValue(1.0,2.0,3.0)
  @test_throws ArgumentError evaluate!(nothing, tpf, v, Jtv)

  # same for k=2
  h = ThirdOrderTensorValue{3,Cov(3),Cov(3),Float64,27}(ntuple(i->Float64(i),27)...)
  @test_throws ArgumentError evaluate!(nothing, tpf, h, Jtv)

  # a fully varianced derivative, as `partial_derivatives` produces, goes
  # through: every index knows how it transforms
  gv = TensorValue{Cov(3),Cov(3),Float64,9}(ntuple(i->Float64(i),9)...)
  got = evaluate!(nothing, tpf, gv, Jtv)
  @test got isa TensorValue{Cov(3),Cov(3)}
  @test collect(evaluate!(nothing, inverse_map(tpf), got, Jtv).data) ≈ collect(gv.data)
end

# ─── BROKEN: differentiating a form-valued FE function ───────────────────────
# `∇` is refused by the testset above, and the `PartialDerivatives` chain rule
# through a `Pushforward` is not implemented yet, so neither spelling works on
# a pushed-forward field. Being addressed separately; see `ReferenceFEs/
# Pullbacks.jl`.

# ─── The AFW moment, pinned by hand ──────────────────────────────────────────
# The lowest order case where the value of a DOF is known in closed form.
# σ_e(dx¹) = ∫_e ι_e^*(dx¹) = x¹(end) − x¹(start). This is the one test that
# detects a constant rescaling of the moment on a single face: the hypotenuse
# of TRI has measure √2, so a double count of the face measure shows up there
# and nowhere else. Neither duality nor trace locality would see it.

@testset "Whitney P⁻₁Λ¹ moments have their exact closed form values" begin
  reffe = pλ(TRI, 1, 1, :P⁻)
  @test num_dofs(reffe) == 3

  dofs = get_dof_basis(reffe)
  vs   = get_vertex_coordinates(TRI)
  e2v  = get_faces(TRI,1,0)

  for (comp, ω) in ((1, ExteriorFormValue{1,Cov(2)}(1.0,0.0)),
                    (2, ExteriorFormValue{1,Cov(2)}(0.0,1.0)))
    got  = evaluate(dofs, ConstantField(ω))
    want = [ vs[e2v[e][2]][comp] - vs[e2v[e][1]][comp] for e in 1:3 ]
    @test got ≈ want
  end
end

# ─── Reference level validation over the whole family ────────────────────────

@testset "Form valued FEEC element: construction, duality and trace locality" begin
  for D in (2,3), r in 1:3, k in 0:D, F in (:P,:P⁻)
    p = simplex(D)
    reffe = pλ(p, r, k, F)
    prebasis = get_prebasis(reffe)

    # the form degree is readable from the element itself
    @test get_form_degree(reffe) == k
    @test Conformity(reffe) isa PΛConformity

    # the moments are unisolvent on the prebasis
    @test num_dofs(reffe) == length(prebasis)

    M = evaluate(get_dof_basis(reffe), get_shapefuns(reffe))
    @test M ≈ Matrix(I, num_dofs(reffe), num_dofs(reffe))

    # Trace locality. An AFW DOF of a face f only sees the trace ι_f^*ω, and a
    # bubble owned by a face g has zero trace on any face not containing g. So
    # σ_f(φ_g) must vanish whenever g ⊄ f. Unlike the duality above -- which is
    # tautological, the shape functions being computed as the dual of the dofs
    # -- this actually exercises the pullback, the pairing and the wiring.
    faces = get_faces(p)
    face_own_dofs = get_face_own_dofs(reffe)
    face_own_funs = get_face_own_funs(prebasis, p, PΛConformity())
    A = evaluate(get_dof_basis(reffe), prebasis)
    for (f, dof_ids) in enumerate(face_own_dofs), (g, fun_ids) in enumerate(face_own_funs)
      (isempty(dof_ids) || isempty(fun_ids) || faces[g] ⊆ faces[f]) && continue
      @test all(iszero, round.(A[dof_ids, fun_ids], digits=10))
    end
  end
end

@testset "test_reference_fe" begin
  for D in (2,3), r in 1:2, k in 0:D, F in (:P,:P⁻)
    test_reference_fe(pλ(simplex(D), r, k, F))
  end
end

# ─── Value type entry point, and parity with the vector proxy elements ───────
# The form valued element is reached by giving an `ExteriorFormValue` type to
# the usual named constructor -- which is what lets a 1-form space be built
# over a mesh of mixed element types, where no single family symbol applies.
#
# Both representations must describe the very same space, so their total and
# per-face DOF counts must agree at equal name and equal `order`. This is what
# validates the AFW periodic table bookkeeping
#   PᵣΛᵏ  on a d-face → P⁻_{r+k−d}   Λ^{d−k}(f)
#   P⁻ᵣΛᵏ on a d-face → P_{r+k−d−1} Λ^{d−k}(f)
# against the existing elements, independently of this implementation.

@testset "Form value type selects the form valued variant of each name" begin
  for D in (2,3)
    p = simplex(D)
    for (name, k, orders) in ((nedelec,        1,   0:2),
                              (nedelec2,       1,   1:3),
                              (raviart_thomas, D-1, 0:2),
                              (bdm,            D-1, 1:3))
      for order in orders
        form  = ReferenceFE(p, name, ExteriorFormValue{k,D,Float64}, order)
        proxy = ReferenceFE(p, name, Float64, order)

        @test get_name(form) == get_name(proxy)
        @test get_form_degree(form) == k
        @test num_dofs(form) == num_dofs(proxy)
        @test map(length, get_face_own_dofs(form)) == map(length, get_face_own_dofs(proxy))
        @test Pushforward(get_name(form), Conformity(form)) isa TensorPushforward
      end
    end

    # The top of the table: D-forms are L² densities whose coefficients span
    # 𝓟_order, and -- unlike their scalar counterpart -- they are transported.
    for order in 0:3
      top = ReferenceFE(p, modal_lagrangian, ExteriorFormValue{D,D,Float64}, order)
      @test num_dofs(top) == binomial(order+D,D)
      @test get_form_degree(top) == D
      @test Conformity(top) isa PΛConformity
    end
  end

  # In 2D the 1-forms and the (D-1)-forms are the same space, so the two names
  # of the trimmed family build the very same element.
  @test num_dofs(ReferenceFE(TRI, nedelec, ExteriorFormValue{1,2,Float64}, 1)) ==
        num_dofs(ReferenceFE(TRI, raviart_thomas, ExteriorFormValue{1,2,Float64}, 1))

  # A name stands for a space, hence for one form degree, and for one dimension
  @test_throws AssertionError ReferenceFE(TRI, nedelec, ExteriorFormValue{2,2,Float64}, 1)
  @test_throws AssertionError ReferenceFE(TET, raviart_thomas, ExteriorFormValue{1,3,Float64}, 1)
  @test_throws AssertionError ReferenceFE(TRI, nedelec, ExteriorFormValue{1,3,Float64}, 1)
  # n-cubes have no form valued polynomial basis for 0 < k < D yet
  @test_throws AssertionError ReferenceFE(QUAD, nedelec, ExteriorFormValue{1,2,Float64}, 1)

  # Only a form value type has a form degree. A `Float64` is not a 0-form --
  # `∧`, `flat_hodge_star` and `pullback` are not defined on it -- and a vector
  # proxy is ambiguous between Λ¹ and Λᴰ⁻¹, which is the whole point.
  @test get_form_degree(ExteriorFormValue{2,3,Float64}) == 2
  @test get_form_degree(ExteriorFormValue{0,3,Float64}) == 0
  @test_throws ErrorException get_form_degree(VectorValue{3,Float64})
  @test_throws ErrorException get_form_degree(Float64)
end

# ─── The FEEC factory hop ────────────────────────────────────────────────────

@testset "ReferenceFE(p,F,r,k,T; DG_calc=true) hops to the right name" begin
  for D in (2,3), (F,k,r,name) in ((:P⁻,1,1,nedelec), (:P,1,1,nedelec2),
                                   (:P⁻,D-1,2,raviart_thomas), (:P,D-1,2,bdm),
                                   (:P⁻,D,2,modal_lagrangian), (:P,D,2,modal_lagrangian))
    p = simplex(D)
    reffe = pλ(p, r, k, F)
    # in 2D, k=1=D-1 and the 1-form branch wins, as `rotate_90` -- the flag
    # that disambiguates the two vector proxies -- has no form valued meaning
    expected = (D == 2 && k == 1) ? (F == :P ? nedelec2 : nedelec) : name
    @test get_name(reffe) == expected
    @test get_form_degree(reffe) == k
    @test Conformity(reffe) isa PΛConformity
  end

  # The form value type may also be given to the factory directly, `k` then
  # being redundant with it.
  for D in (2,3), k in 1:D
    V = ExteriorFormValue{k,D,Float64}
    @test num_dofs(ReferenceFE(simplex(D), :P, 2, k, V)) == num_dofs(pλ(simplex(D), 2, k, :P))
  end

  # A 0-form is NOT a `Float64`: `∧`, `flat_hodge_star` and `pullback` are
  # defined on `ExteriorFormValue{0,D,T}` and not on a real number. So k=0 is a
  # form valued element of its own, with the same DOFs as the scalar one.
  for D in (2,3), F in (:P,:P⁻)
    p = simplex(D)
    reffe = pλ(p, 2, 0, F)
    scalar = ReferenceFE(p, modal_lagrangian, Float64, 2)
    @test get_name(reffe) == modal_lagrangian
    @test Conformity(reffe) isa PΛConformity
    @test get_form_degree(reffe) == 0
    @test num_dofs(reffe) == num_dofs(scalar)
    @test map(length, get_face_own_dofs(reffe)) == map(length, get_face_own_dofs(scalar))
    @test Pushforward(get_name(reffe), Conformity(reffe)) isa TensorPushforward
  end

  # Nodal DOFs are meaningless for `0 < k < D`
  @test_throws AssertionError ReferenceFE(TET, :P, 2, 1, Float64; DG_calc=true, nodal=true)

  # The family still has to match the polytope, at every k -- the form valued
  # branch returns before the proxy table that otherwise checks this
  @test_throws AssertionError ReferenceFE(QUAD, :P, 2, 0, Float64; DG_calc=true, nodal=true)
  @test_throws AssertionError ReferenceFE(TRI, :Q⁻, 2, 0, Float64; DG_calc=true, nodal=true)
  @test_throws AssertionError ReferenceFE(QUAD, :P⁻, 2, 2, Float64; DG_calc=true, nodal=true)

  # On the form valued path a tensor `T` is not a Cartesian product request,
  # and a form type must agree with the `(k,D)` asked for
  @test_throws AssertionError ReferenceFE(TET, :P, 2, 1, VectorValue{3,Float64}; DG_calc=true)
  @test_throws AssertionError ReferenceFE(TET, :P, 2, 1, ExteriorFormValue{2,3,Float64})
  @test_throws AssertionError ReferenceFE(TET, :P, 2, 1, ExteriorFormValue{1,2,Float64})
  # the vector proxy path is unchanged
  @test num_dofs(ReferenceFE(TET, :P, 2, 0, VectorValue{3,Float64})) == 30
end

# ─── The nodal names, on simplices AND n-cubes ───────────────────────────────
# `lagrangian` and `serendipity` cover the two scalar ends of the periodic
# table -- the 𝓟ᵣ⁽⁻⁾Λ⁰/𝓠ᵣ⁻Λ⁰/𝓢ᵣΛ⁰ and 𝓟ᵣ⁽⁻⁾Λᴰ/𝓠ᵣ⁻Λᴰ/𝓢ᵣΛᴰ rows of
# `docs/src/modules/ReferenceFEs.md`. A 0-form and a D-form carry a single
# component, so the nodal element is the scalar one with a form value type; it
# differs only in being tagged `PΛConformity`, hence transported.

@testset "Nodal names carry 0-forms and D-forms, on every polytope" begin
  for (p,D) in ((TRI,2),(TET,3),(QUAD,2),(HEX,3)), order in 1:3
    for k in (0,D)
      V = ExteriorFormValue{k,D,Float64}
      form   = ReferenceFE(p, lagrangian, V, order)
      scalar = ReferenceFE(p, lagrangian, Float64, order)

      @test get_form_degree(form) == k
      @test num_dofs(form) == num_dofs(scalar)
      @test Conformity(form) isa PΛConformity
      @test Pushforward(get_name(form), Conformity(form)) isa TensorPushforward
      test_reference_fe(form)

      if iszero(k)
        # a 0-form has the same trace as a scalar field, hence the same gluing
        @test map(length, get_face_own_dofs(form)) == map(length, get_face_own_dofs(scalar))
      else
        # a D-form has no trace on a facet: every DOF belongs to the cell
        own = get_face_own_dofs(form)
        @test all(isempty, own[1:end-1])
        @test length(own[end]) == num_dofs(form)
        @test own == get_face_own_dofs(form, Conformity(form))
        # The permutations FESpaces consumes: nothing below the cell is glued, so every
        # entry there is empty -- but there must be one *per vertex permutation* of the
        # face, because `CellDofsNonOriented` indexes before testing for emptiness. This
        # is what lets the D-forms live on a mesh whose cells are not sorted.
        perms = get_face_own_dofs_permutations(form)
        @test all(pdofs -> all(isempty, pdofs), perms[1:end-1])
        @test map(length, perms[1:end-1]) ==
              map(length, get_face_vertex_permutations(p)[1:end-1])
        @test only(perms[end]) == collect(1:num_dofs(form))
      end
    end
  end

  # Order 0 at k=D: the piecewise constants, the most used D-form element. This
  # is also where the scalar code takes its `all(orders .== 0) → L2Conformity`
  # branch, which the form valued one must take precedence over.
  for (p,D) in ((TRI,2),(TET,3),(QUAD,2),(HEX,3))
    form = ReferenceFE(p, lagrangian, ExteriorFormValue{D,D,Float64}, 0)
    @test num_dofs(form) == 1
    @test Conformity(form) isa PΛConformity
    @test Pushforward(get_name(form), Conformity(form)) isa TensorPushforward
    @test get_face_own_dofs(form)[end] == [1]
    test_reference_fe(form)
  end

  # Serendipity covers 𝓢ᵣΛ⁰ and 𝓢ᵣΛᴰ, on n-cubes
  for (p,D) in ((QUAD,2),(HEX,3)), order in 2:3
    form = ReferenceFE(p, serendipity, ExteriorFormValue{0,D,Float64}, order)
    @test num_dofs(form) == num_dofs(ReferenceFE(p, serendipity, Float64, order))
    @test Conformity(form) isa PΛConformity
    @test Pushforward(get_name(form), Conformity(form)) isa TensorPushforward
    test_reference_fe(form)

    # 𝓢ᵣΛᴰ ≡ 𝓟ᵣ: the serendipity D-forms are the FULL polynomial densities,
    # not the serendipity ones, so the named constructor must agree with the
    # factory rather than build 𝓢ᵣ
    top = ReferenceFE(p, serendipity, ExteriorFormValue{D,D,Float64}, order)
    @test num_dofs(top) == binomial(order+D,D)
    @test num_dofs(top) == num_dofs(ReferenceFE(p, :S, order, D, Float64; DG_calc=true, nodal=true))
    @test Conformity(top) isa PΛConformity
    test_reference_fe(top)
  end

  # The same elements through the FEEC factory, `nodal=true`
  for (p,D,F) in ((TRI,2,:P),(TET,3,:P⁻),(QUAD,2,:Q⁻),(HEX,3,:Q⁻)), k in (0,D)
    reffe = ReferenceFE(p, F, 2, k, Float64; DG_calc=true, nodal=true)
    @test get_form_degree(reffe) == k
    @test Conformity(reffe) isa PΛConformity
    # 𝓟ᵣ⁻Λᴰ = 𝓟ᵣ₋₁ and 𝓠ᵣ⁻Λᴰ = 𝓠ᵣ₋₁, one order below 𝓟ᵣΛ⁰/𝓠ᵣΛ⁰
    scalar_order = (iszero(k) || F == :P) ? 2 : 1
    @test num_dofs(reffe) == num_dofs(ReferenceFE(p, lagrangian, Float64, scalar_order))
  end

  # 𝓢ᵣΛᴰ ≡ 𝓟ᵣ: the `space=:P` nodal Lagrangian element on an n-cube
  for (p,D) in ((QUAD,2),(HEX,3))
    reffe = ReferenceFE(p, :S, 2, D, Float64; DG_calc=true, nodal=true)
    @test num_dofs(reffe) == binomial(2+D,D)
    @test Conformity(reffe) isa PΛConformity
    @test get_form_degree(reffe) == D
  end
end

# ─── The real pushforward path, through the Field machinery ─────────────────
# Everything above calls `evaluate!` pointwise on hand built tensors. This
# drives the map the way an FE space would: as an `Operation` over the actual
# shape functions. The jacobian is the varianced one that
# `partial_derivatives` of a cell map produces, since that is the only kind
# `TensorPushforward` accepts.
#
# BROKEN: derivatives of the pushed-forward basis are not exercised. The
# `Broadcasting(∇)` reroute of `Pullbacks.jl` feeds the map a half varianced
# reference gradient, which is now refused, and the `PartialDerivatives` chain
# rule through a `Pushforward` is not implemented yet. Being addressed
# separately.

@testset "TensorPushforward over shape function fields" begin
  tpf = TensorPushforward()
  for D in (2,3), k in 0:1, F in (:P,:P⁻)
    reffe = pλ(simplex(D), 2, k, F)
    sf = get_shapefuns(reffe)
    Jt_val = D == 2 ? TensorValue{Cov(2),Con(2)}(2.0,0.3, 0.2,1.5) :
                      TensorValue{Cov(3),Con(3)}(2.0,0.3,0.1, 0.2,1.5,0.4, 0.05,0.2,3.0)
    Jt = ConstantField(Jt_val)
    pts = [ Point(ntuple(_->0.17, D)...), Point(ntuple(i->0.1*i, D)...) ]

    pf = Broadcasting(Operation(tpf))(sf, Jt)

    # values
    want = map(v -> evaluate!(nothing, tpf, v, Jt_val), evaluate(sf, pts))
    @test evaluate(pf, pts) ≈ want

    # the `Broadcasting(∇)` reroute hands the map a reference gradient whose
    # derivation index carries no variance, which `pushforward` refuses
    @test_throws ArgumentError evaluate(Broadcasting(∇)(pf), pts)

    # `∂` is the spelling that must work instead, once the `PartialDerivatives`
    # chain rule through a `Pushforward` lands. Flips to "Unexpected Pass" then,
    # which is the signal to turn it into a real assertion.
    @test_broken evaluate(Broadcasting(∂)(pf), pts) isa AbstractArray

    # DIV must refuse rather than silently produce nonsense through the
    # `LazyArray{<:Fill}` catch-all
    lz = lazy_map(tpf, Fill(sf,1), Fill(Jt,1))
    @test_throws ErrorException Fields.DIV(lz)
  end
end

# ─── DOF scaling ─────────────────────────────────────────────────────────────
# A k-form scales like hᵏ, `k` being read off the value type of the prebasis.
# Were the `TensorPushforward` method not selected, the generic one would reach
# `_scaling_function(::TensorPushforward,D)` and throw, so this also checks
# dispatch.

@testset "get_dofscale_setter_function scales a k-form like h^k" begin
  for D in (2,3), k in 0:D
    reffe = pλ(simplex(D), 2, k, :P)
    setter = get_dofscale_setter_function(reffe, TensorPushforward())
    dofscale = zeros(num_dofs(reffe))
    face_own_dofs = get_face_own_dofs(reffe)
    face_meshsize = fill(2.0, num_faces(get_polytope(reffe)))
    setter(dofscale, face_own_dofs, face_meshsize)
    @test all(≈(2.0^k), dofscale)
  end
end

end # module
