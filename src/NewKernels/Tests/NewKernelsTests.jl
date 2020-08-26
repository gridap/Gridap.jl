module NewKernelsTests

using Test
using Gridap.Arrays: CachedArray
using Gridap.NewKernels
using Gridap.TensorValues
# using LinearAlgebra
# using Gridap.Inference

test_kernel(FunctionKernel(+),(3,2),5)
@test NewKernels.return_types((FunctionKernel(+),FunctionKernel(/)),1,1) == (Int,Float64)

f = BroadcastKernel(+)
a = rand(3,2)
b = 3
c = a .+ b
NewKernels.test_kernel(f,(a,b),c)


k = BroadcastKernel(-)
test_kernel(k,(1,),-1)
test_kernel(k,([1,2],),[-1,-2])
test_kernel(k,(1,2),-1)
test_kernel(k,(1.0,2),-1.0)
test_kernel(k,(1,2.0),-1.0)
test_kernel(k,([1,2],2),[-1,0])
test_kernel(k,(2,[1,2]),[1,0])
test_kernel(k,([3,4],[1,2]),[2,2])

f = BroadcastKernel(⋅)
a = fill(TensorValue(2,0,0,0,2,0,0,0,2),2)
b = VectorValue(1,2,3)
c = zeros(VectorValue{3,Int},2)
broadcast!(⋅,c,a,b)
test_kernel(f,(a,b),c)

end # module
