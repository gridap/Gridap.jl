using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.Fields
using Gridap.FESpaces
using Gridap.Polynomials
using Test

# Start with a PhysicalSpaceCellBasis

# domain = (0,1)
# partition = (3,)
domain = (0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)
order = 1
# order = 2

trian = get_triangulation(model)
quad = CellQuadrature(trian,order)
q = get_coordinates(quad)

polytopes = get_polytopes(model)
cell_to_ctype = get_cell_type(model)

grid = get_grid(model)
cell_map = get_cell_map(grid)

# Test against the ref approach...

# T = VectorValue{2,Float64}
T = Float64

reffes = [LagrangianRefFE(T,p,order) for p in polytopes]

psfs, x  = Gridap.FESpaces.compute_cell_space_physical_space_lagrangian(reffes, cell_to_ctype, cell_map)
sfs, x  = Gridap.FESpaces.compute_cell_space(reffes, cell_to_ctype, cell_map)

r = evaluate(sfs,q)
rg = evaluate(gradient(sfs),q)
rp = evaluate(psfs,q)
rgp = evaluate(gradient(psfs),q)

@test all([ rg[i] ≈ rgp[i] for i in 1:length(rg) ])
@test all([ r[i] ≈ rp[i] for i in 1:length(rg) ])

##
reffes = [RaviartThomasRefFE(T,p,order) for p in polytopes]

psfs, dofp  = Gridap.FESpaces.compute_cell_space_physical_space_moment(reffes, cell_to_ctype, cell_map)
sfs, dof  = Gridap.FESpaces.compute_cell_space(reffes, cell_to_ctype, cell_map)

r = evaluate(sfs,q)
rg = evaluate(gradient(sfs),q)
rp = evaluate(psfs,q)
rgp = evaluate(gradient(psfs),q)

@test all([ r[i] ≈ rp[i] for i in 1:length(rg) ])
@test all([ rg[i] ≈ rgp[i] for i in 1:length(rg) ])

dofp[2]
dof
q[1]
r[1]
rp[1]

##
dof_bases = map(get_dof_basis,reffes)

ctype_to_refnodes = map(get_nodes,dof_bases)
cell_to_refnodes = CompressedArray(ctype_to_refnodes,cell_to_ctype)
cell_physnodes = evaluate(cell_map,cell_to_refnodes)

# Not efficient, create a Kernel
ct_face_moments = map(ReferenceFEs.get_face_moments,dof_bases)
c_face_moments = CompressedArray(ct_face_moments,cell_to_ctype)
ct_face_nodes_dofs = map(ReferenceFEs.get_face_nodes_dofs,dof_bases)
c_face_nodes_dofs = CompressedArray(ct_face_nodes_dofs,cell_to_ctype)
cell_dof_basis = apply( (n,m,nd) -> ReferenceFEs.MomentBasedDofBasis(n,m,nd),
                  cell_physnodes, c_face_moments, c_face_nodes_dofs)

prebasis =  map(get_prebasis,reffes)
cell_prebasis = CompressedArray(prebasis,cell_to_ctype)

cell_matrix = evaluate_dof_array(cell_dof_basis,cell_prebasis)
cell_matrix_inv = apply(inv,cell_matrix)
cell_shapefuns_phys = apply(change_basis,cell_prebasis,cell_matrix_inv)
cell_shapefuns = compose(cell_shapefuns_phys,cell_map)

(cell_shapefuns, cell_dof_basis)
##

##
# If I want new evaluation...
function kernel_evaluate(k::typeof{change_basis},x,cell_prebasis,cell_matrix_inv)
   cell_prebasis_x = evaluate_field_array(cell_prebasis,x)
  apply(mul,cell_prebasis_x,cell_prebasis,cell_matrix_inv)
end
function apply_gradient(k::typeof(change_basis),cell_prebasis,cell_matrix_inv)
   cell_prebasis_grad = gradient(cell_prebasis)
   apply(change_basis,cell_prebasis_grad,cell_matrix_inv)
end
##
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


a1 = Gridap.Arrays.Fill(1.0,3)
b1 = Gridap.Arrays.Fill(1.0,3)
c1 = Gridap.Arrays.Fill(1.0,3)
f(a,b,c) = a+b+c
p1 = Gridap.Arrays.pair_arrays(a1,b1)
p2 = Gridap.Arrays.pair_arrays(p1,c1)

apply(f,a1,b1,c1)
