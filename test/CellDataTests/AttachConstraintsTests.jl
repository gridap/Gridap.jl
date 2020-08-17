module AttachConstraintsTests

using Gridap.Arrays
using Gridap.CellData

ncells = 10
ndofs = 3
ndofs_c = 4
cellmat = [ rand(ndofs,ndofs) for cell in 1:ncells]
cellvec = [ rand(ndofs) for cell in 1:ncells]
cellmatvec = pair_arrays(cellmat,cellvec)
cellconstr = [ rand(ndofs_c,ndofs) for cell in 1:ncells]

a = attach_constraints_rows(cellvec,cellconstr)
r = map( (vec,constr) -> constr*vec ,cellvec,cellconstr)
test_array(a,r)

a = attach_constraints_rows(cellmat,cellconstr)
r = map( (mat,constr) -> constr*mat ,cellmat,cellconstr)
test_array(a,r)

a = attach_constraints_rows(cellmatvec,cellconstr)
r = map( (mat,vec,constr) -> (constr*mat,constr*vec),cellmat,cellvec,cellconstr)
test_array(a,r)

a = attach_constraints_cols(cellmat,cellconstr)
r = map( (mat,constr) -> mat*transpose(constr) ,cellmat,cellconstr)
test_array(a,r)

a = attach_constraints_cols(cellmatvec,cellconstr)
r = map( (mat,vec,constr) -> (mat*transpose(constr),vec),cellmat,cellvec,cellconstr)
test_array(a,r)

#cache = array_cache(a)
#using BenchmarkTools
#@btime getindex!($cache,$a,2)

end # module
