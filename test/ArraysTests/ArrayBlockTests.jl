module ArrayBlockTests

using Test
using Gridap
using Gridap.Arrays
using LinearAlgebra

A = ArrayBlock(reshape([fill(1.0,2,3),fill(0.0,0,0),fill(2.0,2,5), fill(3.0,4,5)],(2,2)), Bool[1 1; 0 1])
B = ArrayBlock([fill(4.0,3),fill(5.0,5)], Bool[1,1])

C = A*B
Ct = transpose(C)
Ct * C

C = A*transpose(A)
mul!(C,A,transpose(A))
mul!(C,A,transpose(A),2.0,1.5)

cache = return_cache(*,A,B)
C = evaluate!(cache,*,A,B)

k = MulAddMap(2.0,1.5)
cache = return_cache(k,A,B,C)
D = evaluate!(cache,k,A,B,C)
@test D == D - D + D

end
