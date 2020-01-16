module SingleFieldAssemblersTests

using Test
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Integration
using Gridap.Fields
using SparseArrays

#using Gridap.FESpaces

include("../../../src/new/FESpaces/FESpaces.jl")
using .FESpaces

domain =(0,1,0,1)
partition = (2,2)
model = CartesianDiscreteModel(domain,partition)

order = 2
grid_topology = get_grid_topology(model)
polytopes = get_polytopes(grid_topology)
reffes = [LagrangianRefFE(Float64,p,order) for p in polytopes]

dirichlet_tags = ["tag_1","tag_6"]
V = GradConformingFESpace(reffes,model,dirichlet_tags)

U = TrialFESpace(V,[4,3])

v = get_cell_basis(V)
u = get_cell_basis(U)

trian = get_triangulation(model)
degree = 2
quad = CellQuadrature(trian,degree)

btrian = BoundaryTriangulation(model)
bquad = CellQuadrature(btrian,degree)

bu = restrict(u,btrian)
bv = restrict(v,btrian)

cellmat = integrate(v*u,trian,quad)
cellvec = integrate(v*3,trian,quad)
cellids = collect(1:num_cells(trian))

bcellmat = integrate(bv*bu,btrian,bquad)
bcellvec = integrate(bv*3,btrian,bquad)
bcellids = get_face_to_cell(btrian)

term_to_cellmat = [cellmat, bcellmat]
term_to_cellvec = [cellvec, bcellvec]
term_to_cellids = [cellids, bcellids]

assem = SparseMatrixAssembler(SparseMatrixCSC,V,U)
#test_assembler(assem,term_to_cellmat,term_to_cellvec,term_to_cellids,term_to_cellids)


end # module
