module TransientFESpacesTests

using Test

using Gridap
using Gridap.Fields
using Gridap.ODEs

u1(t) = x -> (x[1] + x[2]) * t
u2(t) = x -> x[1] * t^2 + x[2] * t

domain = (0, 1, 0, 1)
partition = (5, 5)
model = CartesianDiscreteModel(domain, partition)

order = 1
reffe = ReferenceFE(lagrangian, Float64, order)
V1 = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags="boundary")
U1 = TransientTrialFESpace(V1, u1)
@test test_tfe_space(U1)

V2 = TestFESpace(model, reffe, conformity=:H1, dirichlet_tags=["tag_1", "tag_2"])
U2 = TransientTrialFESpace(V2, [u1, u2])
@test test_tfe_space(U2)

ts = randn(5)
for (U, V, us) in ((U1, V1, [u1]), (U2, V2, [u1, u2]))
  Ut = ∂t(U)
  Utt = ∂tt(U)

  for t in ts
    # Dirichlet values of U
    ust = [ui(t) for ui in us]
    ust = (length(ust) == 1) ? ust[1] : ust
    _U0 = TrialFESpace(V, ust)
    U0 = U(t)

    _ud0 = get_dirichlet_dof_values(_U0)
    ud0 = get_dirichlet_dof_values(U0)

    @test all(ud0 .≈ _ud0)

    # Dirichlet values of ∂t(U)
    ∂tust = [∂t(ui)(t) for ui in us]
    ∂tust = (length(∂tust) == 1) ? ∂tust[1] : ∂tust
    _Ut0 = TrialFESpace(V, ∂tust)
    Ut0 = Ut(t)

    _utd0 = get_dirichlet_dof_values(_Ut0)
    utd0 = get_dirichlet_dof_values(Ut0)

    @test all(utd0 .≈ _utd0)

    # Dirichlet values of ∂tt(U)
    ∂ttust = [∂tt(ui)(t) for ui in us]
    ∂ttust = (length(∂ttust) == 1) ? ∂ttust[1] : ∂ttust
    _Utt0 = TrialFESpace(V, ∂ttust)
    Utt0 = Utt(t)

    _uttd0 = get_dirichlet_dof_values(_Utt0)
    uttd0 = get_dirichlet_dof_values(Utt0)

    @test all(uttd0 .≈ _uttd0)
  end
end

end # module TransientFESpacesTests
