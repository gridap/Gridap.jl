module PΛBasisGradientTests
# Gradient path of RotatingPΛBasis / TrimmedPΛBasis (∇ via Gridap's
# FieldGradientArray machinery) validated against central finite differences
# of the value path. Convention: ∇u[a,j] = ∂u_j/∂x_a (TensorValue{D,D}).

using Gridap.Polynomials
using Gridap.Fields: Point, Broadcasting, ∇
using Gridap.Arrays: evaluate
using Test

test_points(D) = D == 2 ?
    [Point(0.13,0.27), Point(0.3,0.3), Point(0.45,0.45), Point(0.1,0.1)] :
    [Point(0.1,0.2,0.3), Point(0.25,0.25,0.25), Point(0.3,0.1,0.15)]

perturb(x, a, h, D) = Point(ntuple(j -> x[j] + (j == a ? h : 0.0), D)...)

@testset "∇(basis) == finite differences" begin
    h = 1e-5
    for (make, name) in ((RotatingPΛBasis, "untrimmed"), (TrimmedPΛBasis, "trimmed"))
        @testset "$name D=$D r=$r" for (D, rs) in ((2, (1, 2, 3)), (3, (1, 2))), r in rs
            b   = make(Val(D), Float64, r)
            pts = test_points(D)
            vals  = evaluate(b, pts)
            grads = evaluate(Broadcasting(∇)(b), pts)
            for (i, x) in enumerate(pts), w in 1:length(b)
                G = grads[i, w]
                for a in 1:D
                    vp = evaluate(b, [perturb(x, a,  h, D)])[1, w]
                    vm = evaluate(b, [perturb(x, a, -h, D)])[1, w]
                    for j in 1:D
                        fd = (vp.data[j] - vm.data[j]) / (2h)
                        @test isapprox(G[a, j], fd; atol=1e-6)
                    end
                end
            end
        end
    end
end

end # module
