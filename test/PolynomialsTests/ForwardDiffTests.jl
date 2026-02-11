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
b = MonomialBasis{2}(V,order)
evaluate(b,xd)

g = Broadcasting(∇)(b)
evaluate(g,xd)

end # module