module PhysicalBasesTests

using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Arrays
using Gridap.Fields
using Gridap.FESpaces
using Gridap.Polynomials
using Gridap.Integration
using Test

# Start with a PhysicalSpaceCellBasis
a = 1
b = ( a == 1 )

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

T = Float64
reffes = [LagrangianRefFE(T,p,order) for p in polytopes]

dof_bases = map(get_dof_basis,reffes)

FEM = Gridap.FESpaces
cell_dof_basis = FEM._cell_dof_basis_physical_space(dof_bases,cell_to_ctype,cell_map)
cell_dof_basis = GenericCellDofBasis(Val{false}(),cell_dof_basis)

prebasis =  map(get_prebasis,reffes)
cell_prebasis = CompressedArray(prebasis,cell_to_ctype)
cell_prebasis_new = GenericCellBasis(Val{false}(),cell_prebasis,cell_map,Val{false}())

typeof(cell_prebasis)
isa(cell_prebasis,CellBasis)

# cell_matrix = evaluate(cell_dof_basis,cell_prebasis)

# cell_shapefuns = _cell_shape_functions_physical_space(cell_prebasis,cell_dof_basis,cell_map)

psfs, x  = Gridap.FESpaces.compute_cell_space_physical(reffes, cell_to_ctype, cell_map)
sfs, x  = Gridap.FESpaces.compute_cell_space(reffes, cell_to_ctype, cell_map)

# T = VectorValue{2,Float64}
# reffes = [LagrangianRefFE(T,p,order) for p in polytopes]

r = evaluate(sfs,q)
rg = evaluate(gradient(sfs),q)
rp = evaluate(psfs,q)
rgp = evaluate(gradient(psfs),q)

@test all([ rg[i] ≈ rgp[i] for i in 1:length(rg) ])
@test all([ r[i] ≈ rp[i] for i in 1:length(rg) ])

T = VectorValue{2,Float64}
reffes = [LagrangianRefFE(T,p,order) for p in polytopes]

psfs, x  = Gridap.FESpaces.compute_cell_space_physical(reffes, cell_to_ctype, cell_map)
sfs, x  = Gridap.FESpaces.compute_cell_space(reffes, cell_to_ctype, cell_map)

r = evaluate(sfs,q)
rg = evaluate(gradient(sfs),q)
rp = evaluate(psfs,q)
rgp = evaluate(gradient(psfs),q)

@test all([ rg[i] ≈ rgp[i] for i in 1:length(rg) ])
@test all([ r[i] ≈ rp[i] for i in 1:length(rg) ])

func(x) = x

cell_field = convert_to_cell_field(func,cell_map)
isa(cell_field,CellField)

# Now RT elements

T = Float64
order = 0
reffes = [RaviartThomasRefFE(T,p,order) for p in polytopes]

psfs, dofp  = Gridap.FESpaces.compute_cell_space_physical(reffes, cell_to_ctype, cell_map)
sfs, dof  = Gridap.FESpaces.compute_cell_space(reffes, cell_to_ctype, cell_map)

r = evaluate(sfs,q)
rg = evaluate(gradient(sfs),q)
rp = evaluate(psfs,q)
rgp = evaluate(gradient(psfs),q)

@test all([ r[i] ≈ rp[i] for i in 1:length(rp) ])
@test all([ rg[i] ≈ rgp[i] for i in 1:length(rg) ])

end #module

# # If I want new evaluation...
# function kernel_evaluate(k::typeof{change_basis},x,cell_prebasis,cell_matrix_inv)
#    cell_prebasis_x = evaluate_field_array(cell_prebasis,x)
#   apply(mul,cell_prebasis_x,cell_prebasis,cell_matrix_inv)
# end
# function apply_gradient(k::typeof(change_basis),cell_prebasis,cell_matrix_inv)
#    cell_prebasis_grad = gradient(cell_prebasis)
#    apply(change_basis,cell_prebasis_grad,cell_matrix_inv)
# end
# ##
# # Optimisation : evaluate_field_array for AbstractArray with FieldLike
# # Define a new kernel that better treats the inverse
# struct InvKernel <: Kernel end
# function kernel_cache(k::InvKernel,mat)
# end
# function apply_kernel!(cache,k::InvKernel,mat)
# end
# function kernel_cache(k::InvKernel,mat)
# CachedArray(copy(mat))
# end
# function apply_kernel!(cache,k::InvKernel,mat)
#   setsize!(cache,size(mat))
#   m = cache.array
#   fill!(m,zero(m))
#   for i:size(m,1); m[i] = 1; end
#   ldiv!(mat,m)
#   m
# end
# k = InvKernel()
#
# isa(cell_prebasis,CellBasis)
#
# change_basis(cell_prebasis[1],cell_matrix_inv[1])
