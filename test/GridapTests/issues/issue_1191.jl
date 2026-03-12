using Gridap

model = CartesianDiscreteModel((0,1),(2,))
Ω = Interior(model)
mask = [true,false]
Ω1 = Interior(model,mask)
Ω2 = Interior(model,.!mask)

reffe = ReferenceFE(lagrangian,VectorValue{1,Float64},1)
V1 = TestFESpace(Ω1,reffe,dirichlet_tags="boundary")
V2 = TestFESpace(Ω2,reffe,dirichlet_tags="boundary")
U1 = TrialFESpace(V1,VectorValue(0.0))
U2 = TrialFESpace(V2,VectorValue(0.0))

X = MultiFieldFESpace([U1,U2])
Y = MultiFieldFESpace([V1,V2])

u1_basis = get_trial_fe_basis(U1)
u2_basis = get_trial_fe_basis(U2)
v1_basis = get_fe_basis(V1)
v2_basis = get_fe_basis(V2)

op1 = u1_basis⋅v1_basis
op2 = u2_basis⋅v1_basis
op3 = u2_basis⋅v2_basis
op4 = u1_basis⋅v2_basis

@testset "Check all op functions" begin
    ops = [op1, op2, op3, op4]
    for op in ops
        @test try 
            op(Point(0.5,))
            true
        catch
            false
        end
    end
end

u1h = interpolate(x -> VectorValue(2x[1]),U1)
u2h = interpolate(x -> VectorValue(2x[1]),U2)
op1h = u1h⋅v1_basis
op2h = u2h⋅v1_basis
op3h = u2h⋅v2_basis
op4h = u1h⋅v2_basis

@testset "Check all op-h functions" begin
    opsh = [op1h, op2h, op3h, op4h]
    for op in opsh
        @test try 
            op(Point(0.5,))
            true
        catch
            false
        end
    end
end