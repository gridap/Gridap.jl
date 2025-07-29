using Gridap,Gridap.Algebra,Gridap.Arrays,Gridap.Geometry,Gridap.CellData,Gridap.ReferenceFEs,Gridap.FESpaces,Gridap.MultiField

domain = (0,1,0,1)
partition = (5,5)
model = CartesianDiscreteModel(domain,partition)
Γ = BoundaryTriangulation(model)
Λ = SkeletonTriangulation(model)
dΛ = Measure(Λ,2)
V1 = FESpace(Triangulation(model,1:10),ReferenceFE(lagrangian,Float64,1))
V2 = FESpace(model,ReferenceFE(lagrangian,Float64,1))
X = MultiFieldFESpace([V1,V2])
uh = zero(X)

f(xh,yh) = ∫(mean(xh[1])*mean(yh[1])+mean(xh[2])*mean(yh[2]))dΛ
f_jac(xh,dxh,yh) = ∫(mean(dxh[1])*mean(yh[1])+mean(dxh[2])*mean(yh[2]))dΛ

op = FEOperator(f,f_jac,X,X)
jac = jacobian(op,uh)
j = jacobian(uh->f(uh,get_fe_basis(X)),uh)
J = assemble_matrix(j,X,X)

jac == J