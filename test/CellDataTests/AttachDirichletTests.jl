module AttachDirichletTests

using Gridap.Arrays
using Gridap.CellData

ncells = 10
ndofs = 3
cellmat = [ rand(ndofs,ndofs) for cell in 1:ncells]
cellvec = [ rand(ndofs) for cell in 1:ncells]
cellmatvec = pair_arrays(cellmat,cellvec)
cellvals = [ rand(ndofs) for cell in 1:ncells]

cellmatvec_with_dbc = attach_dirichlet(cellmatvec,cellvals)
r = map((mat,vec,vals)-> (mat,vec-mat*vals),cellmat,cellvec,cellvals)
test_array(cellmatvec_with_dbc,r)

#a = cellmatvec_with_dbc
#cache = array_cache(a)
#using BenchmarkTools
#@btime getindex!($cache,$a,2)


end # module
