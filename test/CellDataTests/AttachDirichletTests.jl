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
r = collect(cellmatvec_with_dbc)
test_array(cellmatvec_with_dbc,r)

r = map((mat,vec,vals)-> (mat,vec-mat*vals),cellmat,cellvec,cellvals)
a,b = unpair_arrays(cellmatvec_with_dbc)
ra,rb = unpair_arrays(r) 
test_array(a,ra,≈)
test_array(b,rb,≈)

cellmatvec_with_dbc = attach_dirichlet(cellmat,cellvals)
r = collect(cellmatvec_with_dbc)
test_array(cellmatvec_with_dbc,r)

r = map((mat,vals)-> (mat,-mat*vals),cellmat,cellvals)
a,b = unpair_arrays(cellmatvec_with_dbc)
ra,rb = unpair_arrays(r) 
test_array(a,ra,≈)
test_array(b,rb,≈)

#a = cellmatvec_with_dbc
#cache = array_cache(a)
#using BenchmarkTools
#@btime getindex!($cache,$a,2)


end # module
