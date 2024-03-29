module ZeroMeanFESpacesTests

using Test
using Gridap.Arrays
using Gridap.Geometry
using Gridap.Fields
using Gridap.FESpaces
using Gridap.CellData
using Gridap.ReferenceFEs

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

order = 2

trian = get_triangulation(model)
degree = order
dΩ = Measure(trian,degree)

_V = FESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:L2)

V = ZeroMeanFESpace(_V,dΩ)

cellmat = [rand(4,4) for cell in 1:num_cells(model)]
cellvec = [rand(4) for cell in 1:num_cells(model)]
cellmatvec = pair_arrays(cellmat,cellvec)
test_single_field_fe_space(V,cellmatvec,cellmat,cellvec,trian)

U = TrialFESpace(V)
test_single_field_fe_space(U,cellmatvec,cellmat,cellvec,trian)
@test isa(U,ZeroMeanFESpace)

fun(x) = sin(4*pi*(x[1]+x[2]^2)) + 3
uh = interpolate(fun, U)

mean1 = sum(∫(uh)*dΩ)

tol = 1.0e-10
@test abs(mean1) < tol

V = FESpace(model,ReferenceFE(lagrangian,Float64,order);conformity=:L2,constraint=:zeromean)
@test isa(V,ZeroMeanFESpace)

end # module
