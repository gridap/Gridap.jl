module LawTests

using Test
using Gridap.Arrays
using Gridap.Algebra
using Gridap.TensorValues
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Integration
using Gridap.Fields
using Gridap.FESpaces
using Gridap.CellData

#e = @macroexpand @law function  l(a)
#  println("a is Any") 
#end
#
#display(e)

@law ν(u,x) = (u+1.0)*x[1]
@law dν(du,x) = x[1]*du

domain = (0,1,0,1)
partition = (4,4)
model = CartesianDiscreteModel(domain,partition)

order = 2
diritag = "boundary"
V = GradConformingFESpace(Float64,model,order,diritag)
U = TrialFESpace(V)

trian = get_triangulation(model)
degree = 1
quad = CellQuadrature(trian,degree)
q = get_coordinates(quad)

x = get_physical_coordinate(trian)
u = zero(U)
du = get_cell_basis(U)

r = ν(u,x)
@test isa(r,CellField)
@test evaluate(r,q) == evaluate(operate(ν,u,x),q)

dr = dν(du,x)
@test isa(dr,CellField)
@test is_basis(dr)
@test is_trial(dr)
@test evaluate(dr,q) == evaluate(operate(dν,du,x),q)

end # module
