module MyTests

using Test
using Gridap
using Gridap.Geometry
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.ReferenceFEs
using FillArrays

"""
Keep in mind: There is a function _compose_glues() that 
we could extend! 
"""
sol(x) = sum(x)

D = 2
domain = Tuple(repeat([0,1],D))
cmodel = CartesianDiscreteModel(domain,Tuple(fill(2,D)))
fmodel = refine(cmodel,2)

Ωc = Triangulation(cmodel)
Ωf = Triangulation(fmodel)

Λc = Skeleton(cmodel)
Λf = Skeleton(fmodel)

reffe = ReferenceFE(lagrangian,Float64,1)
Vc = TestFESpace(cmodel,reffe)
Uc = TrialFESpace(Vc)
Vf = TestFESpace(fmodel,reffe)
Uf = TrialFESpace(Vf)

"""
 So this works out of the box for any change between a Dc==Dp triangulation and 
 any other triangulation for which the change is possible. 
"""
v_Ωc = interpolate(sol,Vc)
aux  = Triangulation(ReferenceFE{num_cell_dims(Ωc)},get_adapted_model(Λf))
v_Ωf = change_domain(v_Ωc.plus,aux,ReferenceDomain())
v_Λf = change_domain(v_Ωf,Λf,ReferenceDomain())

v_Ωf = interpolate(sol,Vf)
aux  = Triangulation(ReferenceFE{num_cell_dims(Ωf)},get_background_model(Λc))
v_Ωc = change_domain(v_Ωf.plus,aux,ReferenceDomain())
v_Λc = change_domain(v_Ωc,Λc,ReferenceDomain())


"""
Attempt at an edge adaptivity glue
"""

glue = get_adaptivity_glue(fmodel)
n2o_cells = glue.n2o_faces_map[3]
o2n_cells = glue.o2n_faces_map

ctopo = get_grid_topology(cmodel)
ftopo = get_grid_topology(fmodel)

c2e_new = get_faces(ftopo,2,1)
c2e_old = get_faces(ctopo,2,1)

rr = glue.refinement_rules[1]
rr_grid = Adaptivity.get_ref_grid(rr)

# [n/e][child_id][flid]->clid
const rrule_f_to_c_lid_2D=Vector{Vector{Vector{UInt8}}}(
  [
   [[1,1,3,1], [1,2,1,4], [3,1,3,2], [1,4,2,4]],  # nodes
   [[1,1,3,1], [1,1,1,4], [1,2,3,1], [1,2,1,4]]   # edges
  ])

# [n/e][child_id][flid]->clid_dim
const rrule_f_to_c_dim_2D=Vector{Vector{Vector{UInt8}}}(
  [
   [[0,1,1,2], [1,0,2,1], [1,2,0,1], [2,1,1,0]],  # nodes
   [[1,2,1,2], [1,2,2,1], [2,1,1,2], [2,1,2,1]]   # edges
  ])

nE = num_edges(fmodel)
n2o_edges = -ones(nE)

iO = 1
new_cells = o2n_cells[iO]
new_edges = c2e_new[new_cells]
old_edges = c2e_old[iO]


n2o_edges[]