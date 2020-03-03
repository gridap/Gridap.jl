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

cell_shapefuns = apply(change_basis,cell_prebasis,cell_matrix_inv)
gradient(cell_shapefuns)

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


# Whereas for the standard approach...
gradient(cell_prebasis)


shapefuns =  map(get_shapefuns,reffes)
refshapefuns = CompressedArray(shapefuns,cell_to_ctype)
cell_shapefuns = attachmap(refshapefuns,cell_map)
gradient(cell_shapefuns)























# cell_shapefuns = attachmap(cell_shapefuns,cell_map)
isa(cell_shapefuns,CellBasis)

cell = 1

ns = get_nodes(cell_dof_basis[cell])

pns = evaluate(cell_map[cell],ns)

p_changeofbasis = evaluate(cell_prebasis[cell],pns)

# So if we want to evaluate the basis or its gradient in a set of points, we can
# do the same as `BasisFromChangeOfBasis` but the difference being that the computation
# of the change matrix is a cell_array.

# I think we should create a new `CellBasis` struct that at the get_index level
# performs all these computations, makes use of cache arrays, and provides a
# basis. Before doing this, I want to discuss with Francesc what can we use
# of the pre-existing arrays, etc.


grad_cpb = gradient(cell_prebasis)

evaluate(grad_cpb[cell],pns)

gradient(refprebasis)
