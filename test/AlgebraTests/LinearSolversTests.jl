module LinearSolversTests

using Gridap

using LinearAlgebra
using SparseArrays

ls = LUSolver()

diff1(M) = [ [1.0 zeros(1,M-1)]; diagm(1=>ones(M-1)) - Matrix(1.0I,M,M) ]

sdiff1(M) = sparse(diff1(M))

function Laplacian(Nx, Ny, Lx, Ly)
   dx = Lx / (Nx+1)
   dy = Ly / (Ny+1)
   Dx = sdiff1(Nx) / dx
   Dy = sdiff1(Ny) / dy
   Ax = Dx' * Dx
   Ay = Dy' * Dy
   return kron(sparse(1.0I,Ny,Ny), Ax) + kron(Ay, sparse(1.0I,Nx,Nx))
end

n = 10
A = Laplacian(n,n,1,1)
x = rand(n^2)
b = A*x

test_linear_solver(ls,A,b,x)

end # module LinearSolversTests
