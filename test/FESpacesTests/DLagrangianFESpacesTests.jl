module DLagrangianFESpacesTests

using Test
using Gridap

model = CartesianDiscreteModel(partition=(2,2))

grid = Grid(model)

facelabels = FaceLabels(model)

node_to_label = labels_on_dim(facelabels,0)

tag_to_labels = facelabels.tag_to_labels

diritags = [1,2,5]
dirimasks = [true,true,true]

fespace = DLagrangianFESpace(
  Float64,grid,node_to_label,tag_to_labels,diritags,dirimasks)

r = [[-1, -2, 1, 2], [-3, -4, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12]]

@test collect(fespace.cell_to_dofs) == r

trian = Triangulation(grid)
quad = CellQuadrature(trian,order=2)

bh = FEBasis(fespace)

uh = zero(fespace)

cellmat = integrate(inner(bh,bh),trian,quad)
cellvec = integrate(inner(bh,uh),trian,quad)

ufun(x) = x[1]

nfree = 12
ndiri = 4

test_fe_space(fespace, nfree, ndiri, cellmat, cellvec, ufun)

uh = interpolate(fespace,ufun)
u = CellField(trian,ufun)
e = u - uh

tol = 1.0e-8
@test sum(integrate(inner(e,e),trian,quad)) < tol

ufun1(x) = 1.0
ufun2(x) = 2.0
ufun5(x) = 5.0

dv = interpolate_diri_values(fespace,[ufun1,ufun2,ufun5])

@test dv == [1.0, 5.0, 5.0, 2.0]

# DLagrangianFESpace for vector values

T = VectorValue{2,Float64}

diritags = [1,2,5]
dirimasks = [(true,false),(false,true),(true,true)]

fespace = DLagrangianFESpace(
  T,grid,node_to_label,tag_to_labels,diritags,dirimasks)

bh = FEBasis(fespace)

uh = zero(fespace)

cellmat = integrate(inner(bh,bh),trian,quad)
cellvec = integrate(inner(bh,uh),trian,quad)

ufun(x) = x

nfree = 26
ndiri = 6

test_fe_space(fespace, nfree, ndiri, cellmat, cellvec, ufun)

uh = interpolate(fespace,ufun)
u = CellField(trian,ufun)
e = u - uh

tol = 1.0e-8
@test sum(integrate(inner(e,e),trian,quad)) < tol

ufun1(x) = VectorValue(1.0,-1.0)
ufun2(x) = VectorValue(2.0,-2.0)
ufun5(x) = VectorValue(5.0,-5.0)

dv = interpolate_diri_values(fespace,[ufun1,ufun2,ufun5])

@test dv == [1.0, 5.0, -5.0, 5.0, -5.0, -2.0]

order = 1
fespace = DLagrangianFESpace(T,model,order,diritags,dirimasks)
test_fe_space(fespace, nfree, ndiri, cellmat, cellvec, ufun)

order = 2
T = VectorValue{2,Float64}
diritags = [1,2,5]
dirimasks = [(true,false),(false,true),(true,true)]

fespace = DLagrangianFESpace(T,model,order,diritags,dirimasks)

@test fespace.dirimasks === dirimasks

fespace = DLagrangianFESpace(T,model,order,diritags)

@test fespace.dirimasks == [(true, true), (true, true), (true, true)]

order = 2
T = Float64
diritags = ["physical_tag_1","physical_tag_2","physical_tag_5"]
dirimasks = [false,false,true]
fespace = DLagrangianFESpace(T,model,order,diritags,dirimasks)
@test fespace.dirimasks === dirimasks
@test fespace.diritags == [1,2,5]

fespace = DLagrangianFESpace(T,model,order,diritags)
@test fespace.dirimasks == [(true,), (true,), (true,)]

fespace = DLagrangianFESpace(T,model,order)
@test fespace.diritags == Int[]
@test fespace.dirimasks == Tuple{Bool}[]
@test fespace.ndiridofs == 0

end # module
