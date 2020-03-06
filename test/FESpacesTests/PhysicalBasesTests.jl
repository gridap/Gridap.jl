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
cell_to_ctype = get_cell_type(model)

grid = get_grid(model)
cell_map = get_cell_map(grid)

newsfs, x  = compute_cell_space_physical_space(reffes, cell_to_ctype, cell_map)
r11 = evaluate(newsfs,q)
r22 = evaluate(gradient(newsfs),q)

# If I want new evaluation...
function kernel_evaluate(k::typeof{change_basis},x,cell_prebasis,cell_matrix_inv)
   cell_prebasis_x = evaluate_field_array(cell_prebasis,x)
  apply(mul,cell_prebasis_x,cell_prebasis,cell_matrix_inv)
end
function apply_gradient(k::typeof(change_basis),cell_prebasis,cell_matrix_inv)
   cell_prebasis_grad = gradient(cell_prebasis)
   apply(change_basis,cell_prebasis_grad,cell_matrix_inv)
end
# Optimisation : evaluate_field_array for AbstractArray with FieldLike
# Define a new kernel that better treats the inverse
struct InvKernel <: Kernel end
function kernel_cache(k::InvKernel,mat)
end
function apply_kernel!(cache,k::InvKernel,mat)
end
function kernel_cache(k::InvKernel,mat)
CachedArray(copy(mat))
end
function apply_kernel!(cache,k::InvKernel,mat)
  setsize!(cache,size(mat))
  m = cache.array
  fill!(m,zero(m))
  for i:size(m,1); m[i] = 1; end
  ldiv!(mat,m)
  m
end
k = InvKernel()

isa(cell_prebasis,CellBasis)

change_basis(cell_prebasis[1],cell_matrix_inv[1])
##

# Juno.@enter gradient(cell_prebasis)
# Juno.@enter evaluate(g_cpb,q)

function compute_cell_space_physical_space(reffes, cell_to_ctype, cell_map)

  # Create new dof_basis with nodes in the physical space
  ctype_to_refnodes= map(get_node_coordinates,reffes)
  cell_to_refnodes = CompressedArray(ctype_to_refnodes,cell_to_ctype)
  cell_physnodes = evaluate(cell_map,cell_to_refnodes)

  dof_basis = map(get_dof_basis,reffes)
  # Not efficient, create a Kernel
  cell_dof_basis = apply( nodes -> LagrangianDofBasis(Float64,nodes), cell_physnodes )

  prebasis =  map(get_prebasis,reffes)
  cell_prebasis = CompressedArray(prebasis,cell_to_ctype)

  cell_matrix = evaluate_dof_array(cell_dof_basis,cell_prebasis)
  cell_matrix_inv = apply(inv,cell_matrix)
  cell_shapefuns_phys = apply(change_basis,cell_prebasis,cell_matrix_inv)
  cell_shapefuns = compose(cell_shapefuns_phys,cell_map)

  (cell_shapefuns, cell_dof_basis)
end
