module ForwardDiffTests

using Test
using ForwardDiff
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.Polynomials

dualise(x::Number,N) = ForwardDiff.Dual(x,ForwardDiff.Partials(ntuple(i -> 0.0, N)))
dualise(x::Point,N) = Point(dualise.(x.data,N))

x = [Point(0.,0.),Point(1.,0.)]
xd = dualise.(x,2)

order = 1
V = Float64
b = MonomialBasis(Val(2),V,order)
evaluate(b,xd)

g = Broadcasting(∇)(b)
evaluate(g,xd)

# Reproducer for issue #1286
d0 = ForwardDiff.Dual{ForwardDiff.Tag{Arrays.autodiff_array_gradient}}(1,1) # fails
tv = TensorValue(d0)

r = CachedMatrix(fill(tv,1,2))
i = 1
s = Polynomials.StaticArraysCore.MVector(d0)
k = 1
V = VectorValue{2, Float64}
Polynomials._cartprod_set_derivative!(r,i,s,k,V)

end # module
