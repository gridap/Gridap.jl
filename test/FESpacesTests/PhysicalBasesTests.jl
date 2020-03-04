using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.Fields
using Gridap.FESpaces
using Gridap.Polynomials

# Start with a PhysicalSpaceCellBasis

domain = (0,1)
partition = (2)
model = CartesianDiscreteModel(domain,partition)
order = 1

trian = get_triangulation(model)

quad = CellQuadrature(trian,order*2)
q = get_coordinates(quad)

polytopes = get_polytopes(model)
T = Float64

# In fact, Lagrangian not needed, minor thing
reffes = [LagrangianRefFE(T,p,order) for p in polytopes]

dof_basis = map(get_dof_basis,reffes)

grid_topology = get_grid_topology(model)

cell_to_ctype = get_cell_type(grid_topology)

cell_dof_basis = CompressedArray(dof_basis,cell_to_ctype)

grid = get_grid(model)

cell_map = get_cell_map(grid)

prebasis =  map(get_prebasis,reffes)

cell_to_ctype = get_cell_type(grid_topology)

refprebasis = CompressedArray(prebasis,cell_to_ctype)

cell_prebasis = attachmap(refprebasis,cell_map)

cell_matrix = evaluate_dof_array(cell_dof_basis,cell_prebasis)

cell_matrix_inv = apply(inv,cell_matrix)

isa(cell_prebasis,CellBasis)

change_basis(cell_prebasis[1],cell_matrix_inv[1])
##

evaluate(cell_prebasis,q)
g_cpb = gradient(cell_prebasis)
# Juno.@enter gradient(cell_prebasis)
# Juno.@enter evaluate(g_cpb,q)

shapefuns =  map(get_shapefuns,reffes)
refshapefuns = CompressedArray(shapefuns,cell_to_ctype)
cell_shapefuns = attachmap(refshapefuns,cell_map)
g_csf = gradient(cell_shapefuns)
# Juno.@enter gradient(cell_shapefuns)
# Juno.@enter evaluate(g_csf,q)
typeof(g_csf)
# Juno.@enter evaluate(g_csf,q)

cell_physshapefuns = apply(change_basis,cell_prebasis,cell_matrix_inv)
# I don't think I need this step, just to check whether it solved the problem below
cell_physshapefuns = attachmap(cell_physshapefuns,cell_map)
g_cpsf = gradient(cell_physshapefuns)
# Juno.@enter evaluate(g_cpsf,q)
evaluate(g_cpsf,q)

# # # k = g_cpsf.g.value
# # # f = g_cpsf.f
# # # g = apply(k,f...)
# # # h = testitem(g)
# # # b = h.field
# # # x = testitems(q...)
# # # Gridap.Polynomials.gradient_cache(b,x)
# #
# # bs1 = g_csf
# #
# # k1 = bs1.g.value
# # f1 = bs1.f
# # g1 = apply(k1,f1...)
# # h1 = testitem(g1)
# # b1 = h1.field
# # x1 = testitems(q...)
# # Gridap.Polynomials.gradient_cache(b1,x1)
#
# bs2 = g_cpsf
#
# k2 = bs2.g.value
# f2 = bs2.f
# g2 = apply(k2,f2...)
# h2 = testitem(g2)
# b2 = h2.field
# x2 = testitems(q...)
# Gridap.Polynomials.gradient_cache(b2,x2)
##

function compute_cell_space_physical_space(reffes, cell_to_ctype, cell_map)

  dof_basis = map(get_dof_basis,reffes)
  cell_dof_basis = CompressedArray(dof_basis,cell_to_ctype)

  prebasis =  map(get_prebasis,reffes)
  refprebasis = CompressedArray(prebasis,cell_to_ctype)
  cell_prebasis = attachmap(refprebasis,cell_map)

  cell_matrix = evaluate_dof_array(cell_dof_basis,cell_prebasis)
  cell_matrix_inv = apply(inv,cell_matrix)
  cell_shapefuns = apply(change_basis,cell_prebasis,cell_matrix_inv)

  # shapefuns =  map(get_shapefuns,reffes)
  # refshapefuns = CompressedArray(shapefuns,cell_to_ctype)
  # cell_shapefuns = attachmap(cell_shapefuns,cell_map)

  (cell_shapefuns, cell_dof_basis)
end
