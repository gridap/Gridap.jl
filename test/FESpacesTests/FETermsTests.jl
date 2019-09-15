module FETermsTests

using Test
using Gridap
using Gridap.SkeletonCellFields: SkeletonCellMatrix
using Gridap.SkeletonCellFields: SkeletonCellVector

ufun(x) = x[1] + x[2]
bfun(x) = x[1]
gfun(x) = x[2]
sfun(x) = x[1]

model = CartesianDiscreteModel(partition=(4,4))

order = 1
diritag = "boundary"
fespace = H1ConformingFESpace(Float64,model,order,diritag)

V = TestFESpace(fespace)
U = TrialFESpace(fespace,ufun)

trian = Triangulation(model)
quad = CellQuadrature(trian,order=2)

tags = [7,6]
btrian = BoundaryTriangulation(model,tags)
bquad = CellQuadrature(btrian,order=2)

tags = [9,]
strian = SkeletonTriangulation(model,tags)
squad = CellQuadrature(strian,order=2)

v = FEBasis(V)
du = FEBasis(U)
uhd = zero(U)
uh = interpolate(U,ufun)

bfield = CellField(trian,bfun)
gfield = CellField(btrian,gfun)
sfield = CellField(strian,sfun)

a(v,u) = inner(v,u)
l(v) = inner(v,bfield)
g(v) = inner(v,gfield)


t1 = AffineFETerm(a,l,trian,quad)

cm = setup_cell_matrix(t1,v,du)
@test isa(cm,CellMatrix)
@test length(cm) == ncells(trian)

cm = setup_cell_jacobian(t1,uh,v,du)
@test isa(cm,CellMatrix)
@test length(cm) == ncells(trian)

cv = setup_cell_vector(t1,v,uhd)
@test isa(cv,CellVector)
@test length(cv) == ncells(trian)

cv = setup_cell_residual(t1,uh,v)
@test isa(cv,CellVector)
@test length(cv) == ncells(trian)

cn = setup_cell_ids(t1)
@test isa(cn,CellNumber)
@test length(cn) == ncells(trian)

t2 = AffineFETerm(a,g,btrian,bquad)

cm = setup_cell_matrix(t2,v,du)
@test isa(cm,CellMatrix)
@test length(cm) == ncells(btrian)

cv = setup_cell_vector(t2,v,uhd)
@test isa(cv,CellVector)
@test length(cv) == ncells(btrian)

cn = setup_cell_ids(t2)
@test isa(cn,CellNumber)
@test length(cn) == ncells(btrian)

sa(v,u) = inner(jump(v),jump(u))
s(v) = inner(jump(v),sfield)

t6 = AffineFETerm(sa,s,strian,squad)

cm = setup_cell_matrix(t6,v,du)
@test isa(cm,SkeletonCellMatrix)

cv = setup_cell_vector(t6,v,uhd)
@test isa(cv,SkeletonCellVector)

cm = setup_cell_jacobian(t6,uh,v,du)
@test isa(cm,SkeletonCellMatrix)

cv = setup_cell_residual(t6,uh,v)
@test isa(cv,SkeletonCellVector)

cn = setup_cell_ids(t6)
@test isa(cn,Tuple{<:CellNumber,<:CellNumber})

t3 = FESource(g,btrian,bquad)

cv = setup_cell_vector(t3,v,uhd)
@test isa(cv,CellVector)
@test length(cv) == ncells(btrian)

cv = setup_cell_residual(t3,uh,v)
@test isa(cv,CellVector)
@test length(cv) == ncells(btrian)

cn = setup_cell_ids(t3)
@test isa(cn,CellNumber)
@test length(cn) == ncells(btrian)

cm = setup_cell_matrix(t3,v,du)
@test cm === nothing

t4 = LinearFETerm(a,trian,quad)

cm = setup_cell_matrix(t4,v,du)
@test isa(cm,CellMatrix)
@test length(cm) == ncells(trian)

cv = setup_cell_vector(t4,v,uhd)
@test isa(cv,CellVector)
@test length(cv) == ncells(trian)

cn = setup_cell_ids(t4)
@test isa(cn,CellNumber)
@test length(cn) == ncells(trian)

jac(uh,v,du) = a(v,du)
res(uh,v) = a(v,uh) - l(v)

t5 = NonLinearFETerm(res,jac,trian,quad)

cm = setup_cell_jacobian(t5,uh,v,du)
@test isa(cm,CellMatrix)
@test length(cm) == ncells(trian)

cv = setup_cell_residual(t5,uh,v)
@test isa(cv,CellVector)
@test length(cv) == ncells(trian)

cn = setup_cell_ids(t5)
@test isa(cn,CellNumber)
@test length(cn) == ncells(trian)

assem = SparseMatrixAssembler(V,U)

cms = setup_cell_jacobian(uh,v,du,t1)

cms = setup_cell_jacobian(uh,v,du,t1,t2,t6)

mat = assemble(assem,cms...)

cvs = setup_cell_residual(uh,v,t1,t2,t6)

vec = assemble(assem,cvs...)

cvs = setup_cell_vector(v,uhd,t1,t2,t6)

vec = assemble(assem,cvs...)

cms = setup_cell_matrix(v,du,t1,t2,t6)

mat = assemble(assem,cms...)

end #module
