module KernelsTests

using Gridap
using TensorValues

k = NumberKernelFromFunction(-)

test_number_kernel(k,1,4,3)

test_number_kernel(k,-1.0,1.0)

k = NumberKernelFromFunction(sum)

test_number_kernel(k,7,[1,2,4])

k = ArrayKernelFromBroadcastedFunction(*)

test_array_kernel(k,[3,4,3],[1,2,3],[3,2,1])

test_array_kernel(k,[6,4,2],2,[3,2,1])

test_array_kernel(k,[6,4,2],[3,2,1],2)

test_array_kernel(k,[6,4,2],[3,2,1],2,1)

u = rand(Int,2,3,1)
v = rand(Int,1,3,4)
w = u .* v
test_array_kernel(k,w,u,v)
test_array_kernel(k,w,u,v,1)

k = ArrayKernelFromBroadcastedFunction(+)

v1 = VectorValue(2,3)
v2 = VectorValue(3,2)
v3 = VectorValue(1,2)
u = [v1,v2,v3]
v = v1
w = u .+ v
test_array_kernel(k,w,u,v)

w = broadcast(+,u,v,0)
test_array_kernel(k,w,u,v,0)


end # module
