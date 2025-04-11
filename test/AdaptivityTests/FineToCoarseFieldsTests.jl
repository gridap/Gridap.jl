module FineToCoarseFieldsTests

using Test
using Gridap
using Gridap.Arrays
using Gridap.Algebra
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using Gridap.FESpaces
using FillArrays

sol(x) = x[1] + x[2]

D = 2
order  = 1
qorder = order*2+1
domain = Tuple(repeat([0,1],D))

parent = CartesianDiscreteModel(domain,Tuple(fill(2,D)))
model  = refine(parent)

trian  = Triangulation(model)
ctrian = Triangulation(parent)

qorder = order*2+1
dΩ_c = Measure(ctrian,qorder)
dΩ_f = Measure(trian,qorder)
dΩ_cf = Measure(ctrian,trian,qorder)

glue = get_adaptivity_glue(model)
rrules = Adaptivity.get_old_cell_refinement_rules(glue)

cell_reffe_lag = lazy_map(rr->ReferenceFE(get_polytope(rr),rr,lagrangian,Float64,order),rrules)
lazy_map(test_reference_fe,cell_reffe_lag)
cell_reffe_ned = lazy_map(rr->ReferenceFE(get_polytope(rr),rr,nedelec,Float64,order),rrules)
lazy_map(test_reference_fe,cell_reffe_ned)

# Lagrangian tests
reffe  = ReferenceFE(lagrangian,Float64,order)
V_c    = TestFESpace(parent,rrules,reffe;conformity=:H1,dirichlet_tags="boundary")
U_c    = TrialFESpace(V_c,sol)

V_f    = TestFESpace(model,reffe;conformity=:H1,dirichlet_tags="boundary")
U_f    = TrialFESpace(V_f,sol)

test_fe_space(U_c)
test_fe_space(U_f)

u_c = interpolate_everywhere(sol,U_c)
u_f = interpolate_everywhere(sol,U_f)

u_fc = interpolate(u_f,U_c)
u_fc2 = interpolate_everywhere(u_f,U_c)

eh = u_c - u_f
@test sum(∫(eh⋅eh)*dΩ_f) < 1.e-12

eh2 = u_c - u_fc
@test sum(∫(eh2⋅eh2)*dΩ_c) < 1.e-12

eh3 = u_c - u_fc2 
@test sum(∫(eh3⋅eh3)*dΩ_c) < 1.e-12

# Moment-based (Nedelec) tests
sol((x,y)) = 2*VectorValue(-y,x)

reffe  = ReferenceFE(nedelec,Float64,order)
V_c    = TestFESpace(parent,rrules,reffe;dirichlet_tags="boundary")
U_c    = TrialFESpace(V_c,sol)

V_f    = TestFESpace(model,reffe;dirichlet_tags="boundary")
U_f    = TrialFESpace(V_f,sol)

test_fe_space(U_c)
test_fe_space(U_f)

u_c = interpolate_everywhere(sol,U_c)
u_f = interpolate_everywhere(sol,U_f)

u_fc = interpolate(u_f,U_c)
u_fc2 = interpolate_everywhere(u_f,U_c)

eh = u_c - u_f
@test sum(∫(eh⋅eh)*dΩ_f) < 1.e-12

eh2 = u_c - u_fc
@test sum(∫(eh2⋅eh2)*dΩ_c) < 1.e-12

eh3 = u_c - u_fc2 
@test sum(∫(eh3⋅eh3)*dΩ_c) < 1.e-12

modelH = CartesianDiscreteModel((0,1,0,1),(1,1))
modelh = refine(modelH,2)
reffe = LagrangianRefFE(Float64,QUAD,1)
XH = TestFESpace(modelH,reffe)
xH = get_fe_basis(XH)
xHh = change_domain(xH,get_triangulation(modelh),ReferenceDomain())
evaluate(Gridap.CellData.get_data(xHh)[1],[Point(0.0,0.0),Point(0.5,0.5)])
evaluate(Gridap.CellData.get_data(xHh)[1],Point(0.5,0.5))

function test_PR_1074()

    function setup_coarse_discrete_model()
        ptr  = [ 1, 5 ]
        data = [ 1,2,3,4  ]
        cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
        node_coordinates = Vector{Point{2,Float64}}(undef,4)
        node_coordinates[1]=Point{2,Float64}(0.0,0.0)
        node_coordinates[2]=Point{2,Float64}(1.0,0.0)
        node_coordinates[3]=Point{2,Float64}(0.0,1.0)
        node_coordinates[4]=Point{2,Float64}(1.0,1.0)

        polytope=QUAD
        scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
        cell_types=collect(Fill(1,length(cell_vertex_lids)))
        cell_reffes=[scalar_reffe]
        grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                                cell_vertex_lids,
                                                cell_reffes,
                                                cell_types,
                                                Gridap.Geometry.NonOriented())
        Gridap.Geometry.UnstructuredDiscreteModel(grid)
    end

    function setup_adapted_discrete_model(parent)
        ptr  = [ 1, 5, 9 ]
        data = [ 1,2,3,4, 2,5,4,6 ] 
        cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
        node_coordinates = Vector{Point{2,Float64}}(undef,6)
        node_coordinates[1]=Point{2,Float64}(0.0,0.5)
        node_coordinates[2]=Point{2,Float64}(0.5,0.5)
        node_coordinates[3]=Point{2,Float64}(0.0,1.0)
        node_coordinates[4]=Point{2,Float64}(0.5,1.0)
        node_coordinates[5]=Point{2,Float64}(1.0,0.5)
        node_coordinates[6]=Point{2,Float64}(1.0,1.0)

        polytope=QUAD
        scalar_reffe=Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
        cell_types=collect(Fill(1,length(cell_vertex_lids)))
        cell_reffes=[scalar_reffe]
        grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                                cell_vertex_lids,
                                                cell_reffes,
                                                cell_types,
                                                Gridap.Geometry.NonOriented())
        model=Gridap.Geometry.UnstructuredDiscreteModel(grid)

        n2o_faces_map=Vector{Vector{Int64}}(undef,3)
        n2o_faces_map[3]=[1,1]
        n2o_cell_to_child_id=[2,3]
        
        ref_rules = [Gridap.Adaptivity.RefinementRule(Gridap.ReferenceFEs.LagrangianRefFE(Float64,QUAD,1),2)]

        glue=Gridap.Adaptivity.AdaptivityGlue(Gridap.Adaptivity.RefinementGlue(),
                                        n2o_faces_map,
                                        n2o_cell_to_child_id,
                                        ref_rules)
        Gridap.Adaptivity.AdaptedDiscreteModel(model,parent,glue)
    end

    coarse_model=setup_coarse_discrete_model()
    fine_model=setup_adapted_discrete_model(coarse_model)

    order=0

    VH = FESpace(coarse_model,
                ReferenceFE(raviart_thomas,Float64,order),
                conformity=:Hdiv)

    UH = TrialFESpace(VH)
        
    Vh = FESpace(fine_model,
                ReferenceFE(raviart_thomas,Float64,order),
                conformity=:Hdiv)
        
    Uh = TrialFESpace(Vh)

    u(x) = VectorValue(x[1],x[2])
    uh  = interpolate(u,Uh)
    uH  = interpolate(u,UH)
    uhH = interpolate(uh,UH) 

    @test get_free_dof_values(uhH)[2] ≈ 1.0
end

test_PR_1074()

end