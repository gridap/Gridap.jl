module ComplexChangeDomainTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using FillArrays

dot_prod(u,v,Ω::Triangulation) = sum(∫(v⋅u)*Measure(Ω,qorder))
dot_prod(u,v,dΩ::Measure) = sum(∫(v⋅u)*dΩ)
norm2(u,Ω::Triangulation) = sum(∫(u⋅u)*Measure(Ω,qorder))
norm2(u,dΩ::Measure) = sum(∫(u⋅u)*dΩ)

order = 1
qorder = 2*order+1
sol(x) = sum(x)

D = 2
domain = Tuple(repeat([0,1],D))
cmodel = CartesianDiscreteModel(domain,Tuple(fill(2,D)))
fmodel = refine(cmodel,2)

Ωc = Triangulation(cmodel)
Ωf = Triangulation(fmodel)

reffe = ReferenceFE(lagrangian,Float64,order)
Vc = TestFESpace(cmodel,reffe)
Uc = TrialFESpace(Vc)
Vf = TestFESpace(fmodel,reffe)
Uf = TrialFESpace(Vf)

"""
  BodyFittedTriangulation Views
"""

Ωc_view = view(Ωc,[1,2])
Ωf_view = view(Ωf,[1,2,3,4,5,6,7,8])

v_Ωf = interpolate(sol,Vf)
v_Ωf_view = change_domain(v_Ωf,Ωf_view,ReferenceDomain())

v_Ωc_view_view = change_domain(v_Ωf_view,Ωc_view,ReferenceDomain())

v_Ωc = interpolate(sol,Vc)
v_Ωc_view = change_domain(v_Ωc,Ωc_view,ReferenceDomain())

v_Ωf_view_view = change_domain(v_Ωc_view,Ωf_view,ReferenceDomain())

@test norm2(v_Ωc_view_view,Ωc_view) ≈ norm2(v_Ωc_view,Ωc_view)
@test norm2(v_Ωf_view_view,Ωf_view) ≈ norm2(v_Ωf_view,Ωf_view)

# Same but automatic
@test norm2(v_Ωc_view,Ωf_view) ≈ norm2(v_Ωf_view,Ωf_view)
@test norm2(v_Ωf_view,Ωc_view) ≈ norm2(v_Ωc_view,Ωc_view)

# Integrate over Ωc_view by integrating first in Ωf and then moving contributions
dΩ_fc = Measure(Ωc_view,Ωf,qorder)
@test dot_prod(v_Ωf,v_Ωc,dΩ_fc) ≈ dot_prod(v_Ωf,v_Ωc,Ωc_view)

"""
  BodyFitted -> Boundary
"""
Γc = Boundary(cmodel)
Γf = Boundary(fmodel)

v_Ωc = interpolate(sol,Vc)
v_Ωf = interpolate(sol,Vf)

v_Γf = change_domain(v_Ωf,Γf,ReferenceDomain())
v_Γc = change_domain(v_Ωc,Γc,ReferenceDomain())

@test norm2(v_Ωf,Γc) ≈ norm2(v_Γc,Γc)
@test norm2(v_Ωc,Γf) ≈ norm2(v_Γf,Γf)

# Integrate in Γf then add contributions for each coarse cell
dΩ_fc = Measure(Ωc,Γf,qorder)
@test dot_prod(v_Ωf,v_Ωc,dΩ_fc) ≈ dot_prod(v_Ωf,v_Ωc,Γf)

"""
  BodyFitted -> Skeleton

  Not currently working, we need to treat SkeletonPairs separately....
"""

Λc = Skeleton(cmodel)
Λf = Skeleton(fmodel)

for sign in [:plus,:minus]
  v_Ωc = getproperty(interpolate(sol,Vc),sign)
  v_Ωf = getproperty(interpolate(sol,Vf),sign)

  v_Λf = change_domain(v_Ωf,Λf,ReferenceDomain())
  v_Λc = change_domain(v_Ωc,Λc,ReferenceDomain())

  @test norm2(v_Ωf,Λc) ≈ norm2(v_Λc,Λc)
  @test norm2(v_Ωc,Λf) ≈ norm2(v_Λf,Λf)
end

# Test that reproduces a BUG spotted from GridapP4est.jl in the following
# issue: https://github.com/gridap/GridapP4est.jl/issues/86
# The code of the test reproduces the error without relying on GridapP4est data 
# structures.

function get_model1()
  m1s_grid_node_coordinates = 
    VectorValue{2, Float64}[(0.0, 0.0), (0.5, 0.0), (0.0, 1.0), (0.5, 1.0), (1.0, 0.0), (1.0, 1.0)]
  m1s_grid_cell_node_ids = 
    Gridap.Arrays.Table(Vector{Int32}[[1, 2, 3, 4], [2, 5, 4, 6]])

  Gridap.Geometry.UnstructuredDiscreteModel(Gridap.Geometry.UnstructuredGrid(
    m1s_grid_node_coordinates,
    m1s_grid_cell_node_ids,
    [Gridap.ReferenceFEs.LagrangianRefFE(Float64,QUAD,1)],
    collect(Fill(Int8(1),2)),
    Gridap.Geometry.Oriented(),
    has_affine_map=true))
end

function get_ref_rules()
  ref_rule_1=Gridap.Adaptivity.RefinementRule(Gridap.Adaptivity.WithoutRefinement(), 
                                    QUAD, 
                                 Gridap.Adaptivity.compute_reference_grid(QUAD,1))

  ref_rule_2=Gridap.Adaptivity.RefinementRule(Gridap.Adaptivity.GenericRefinement(), 
                                              QUAD, 
                                              Gridap.Adaptivity.compute_reference_grid(QUAD,2))
  ref_rule_1, ref_rule_2
end

function get_model2(model1)
    m2s_grid_node_coordinates = 
          VectorValue{2, Float64}[(0.0, 0.0), 
                                  (0.5, 0.0), 
                                  (0.0, 1.0), 
                                  (0.5, 1.0), 
                                  (0.75, 0.0), 
                                  (0.75, 0.5), 
                                  (1.0, 0.0), 
                                  (1.0, 0.5), 
                                  (0.75, 1.0), 
                                  (1.0, 1.0),
                                  (0.5, 0.5)]
    m2s_grid_cell_node_ids = 
        Gridap.Arrays.Table(Vector{Int32}[[1, 2, 3, 4], 
                                          [2, 5, 11, 6], 
                                          [5, 7, 6, 8], 
                                          [11, 6, 4, 9], 
                                          [6, 8, 9, 10]])
     
    ref_rule_1, ref_rule_2 = get_ref_rules()

    ptrs=[1,2]
    values=Vector{Gridap.Adaptivity.RefinementRule}(undef, 2)
    values[1]=ref_rule_1
    values[2]=ref_rule_2
    refinement_rules=Gridap.Arrays.CompressedArray(values, ptrs)

    n2o_faces_map = Vector{Vector{Int}}(undef, 3)
    n2o_faces_map[1] = Vector{Int}(undef, 11)
    n2o_faces_map[2] = Vector{Int}(undef, 16)
    n2o_faces_map[3] = [1, 2, 2, 2, 2]
    n2o_cell_to_child_id=[1, 1, 2, 3, 4]

    glue_m1_to_m2=Gridap.Adaptivity.AdaptivityGlue(Gridap.Adaptivity.RefinementGlue(),
                   n2o_faces_map,
                   n2o_cell_to_child_id,
                   refinement_rules)

    model2=Gridap.Geometry.UnstructuredDiscreteModel(Gridap.Geometry.UnstructuredGrid(
                                                            m2s_grid_node_coordinates,
                                                            m2s_grid_cell_node_ids,
                                                            [Gridap.ReferenceFEs.LagrangianRefFE(Float64,QUAD,1)],
                                                            collect(Fill(Int8(1),5)),
                                                            Gridap.Geometry.Oriented(),
                                                            has_affine_map=true))

    Gridap.Adaptivity.AdaptedDiscreteModel(model2, model1, glue_m1_to_m2)
end 

function get_model3(model2)
    m3s_grid_node_coordinates = 
       VectorValue{2, Float64}[(0.0, 0.0), (0.5, 0.0), (0.0, 1.0), (0.5, 1.0), (1.0, 0.0), (1.0, 1.0)]
    m3s_grid_cell_node_ids = 
        Gridap.Arrays.Table(Vector{Int32}[[1, 2, 3, 4], [2, 5, 4, 6]])

    model3=Gridap.Geometry.UnstructuredDiscreteModel(Gridap.Geometry.UnstructuredGrid(
        m3s_grid_node_coordinates,
        m3s_grid_cell_node_ids,
        [Gridap.ReferenceFEs.LagrangianRefFE(Float64,QUAD,1)],
        collect(Fill(Int8(1),2)),
        Gridap.Geometry.Oriented(),
        has_affine_map=true))

    n2o_faces_map = Vector{Gridap.Arrays.Table{Int64, Vector{Int64}, Vector{Int64}}}(undef, 3)
    n2o_faces_map[3] = Gridap.Arrays.Table(Vector{Int64}[[1], [2, 3, 4, 5]])
    n2o_cell_to_child_id = Gridap.Arrays.Table(Vector{Int64}[[1], [1, 2, 3, 4]])


    ptrs=[1, 2, 2, 2, 2]
    values=Vector{Gridap.Adaptivity.RefinementRule}(undef, 2)
    ref_rule_1, ref_rule_2 = get_ref_rules()
    values[1]=ref_rule_1
    values[2]=ref_rule_2
    refinement_rules=Gridap.Arrays.CompressedArray(values, ptrs)

    glue_m2_to_m3=Gridap.Adaptivity.AdaptivityGlue(Gridap.Adaptivity.MixedGlue(),
                   n2o_faces_map,
                   n2o_cell_to_child_id,
                   refinement_rules)

    Gridap.Adaptivity.AdaptedDiscreteModel(model3, model2, glue_m2_to_m3)
end

model1=get_model1()
model2=get_model2(model1)
model3=get_model3(model2)

t2=Triangulation(ReferenceFE{2}, model2, [false,true,true,true,true])
t3=Triangulation(ReferenceFE{2}, model3,[false,true])

u(x) = x[1]+x[2]

order=1
Vh1 = FESpace(t2,
            ReferenceFE(lagrangian,Float64,order),
            conformity=:H1)

Uh1 = TrialFESpace(Vh1)
uh1 = interpolate(u,Uh1)    

Vh2 = FESpace(t3,
            ReferenceFE(lagrangian,Float64,order),
            conformity=:H1)
    
Uh2 = TrialFESpace(Vh2)

uh2 = interpolate(uh1, Uh2)

dΩ = Measure(t3, 2)
@test sum(∫((uh2-u)⋅(uh2-u))dΩ) < 1e-10

end