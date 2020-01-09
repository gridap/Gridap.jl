module SparseMatricesTests

using Test
using SparseArrays
using Gridap.Algebra

M = SparseMatrixCSC

I = Int[]
J = Int[]
V = Float64[]

push_coo!(M,I,J,V,1,1,3.0)
push_coo!(M,I,J,V,2,1,4.0)
push_coo!(M,I,J,V,2,2,5.0)
push_coo!(M,I,J,V,2,1,6.0)
finalize_coo!(M,I,J,V,2,2)

A = sparse_from_coo(M,I,J,V,2,2)
@test isa(A,M)
@test collect(A) == [3.0 0.0; 10.0 5.0]

end # module
