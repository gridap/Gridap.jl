module NewKernelsTests

using Test
using Gridap.Arrays
using Gridap.NewKernels
using Gridap.TensorValues
using LinearAlgebra

NewKernels.test_kernel(FunctionKernel(+),(3,2),5)
@test NewKernels.return_types((FunctionKernel(+),FunctionKernel(/)),1,1) == (Int,Float64)

f = BroadcastKernel(+)
a = rand(3,2)
b = 3
c = a .+ b
test_kernel(f,(a,b),c)
@test ==(z,y)
@test typeof(z) == return_type(f,x...)
isa(FunctionKernel,NewKernel)

@test kernel_return_types((+,/),1,1) == (Int,Float64)

f = bcast(+)
a = rand(3,2)
b = 3
c = a .+ b
test_kernel(f,(a,b),c)

#test_kernel(1,(),1)
#test_kernel(1,(1,),1)
#
#test_kernel([1,2],(),[1,2])
#test_kernel([1,2],(1,),[1,2])

k = elem(-)
test_kernel(k,(1,),-1)
test_kernel(k,([1,2],),[-1,-2])
test_kernel(k,(1,2),-1)
test_kernel(k,(1.0,2),-1.0)
test_kernel(k,(1,2.0),-1.0)
test_kernel(k,([1,2],2),[-1,0])
test_kernel(k,(2,[1,2]),[1,0])
test_kernel(k,([3,4],[1,2]),[2,2])

k = contract(-)
test_kernel(k,([3,4],[1,2]),3-1+4-2)

f = bcast(⋅)
a = fill(TensorValue(2,0,0,0,2,0,0,0,2),2)
b = VectorValue(1,2,3)
c = zeros(VectorValue{3,Int},2)
broadcast!(⋅,c,a,b)
test_kernel(f,(a,b),c)

end # module
