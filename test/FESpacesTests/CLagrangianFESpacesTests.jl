module CLagrangianFESpacesTests

using Test
using Gridap
using Gridap.CLagrangianFESpaces: grid_from_model_and_order

# CLagrangianFESpace for scalar values

model = CartesianDiscreteModel(partition=(2,2))

grid = Grid(model)

facelabels = FaceLabels(model)

node_to_label = labels_on_dim(facelabels,0)

tag_to_labels = facelabels.tag_to_labels

diritags = [1,2,5]
dirimasks = [true,true,true]

fespace = CLagrangianFESpace(
  Float64,grid,node_to_label,tag_to_labels,diritags,dirimasks)

r = [[-1, -2, 1, 2], [-2, -3, 2, 3], [1, 2, 4, 5], [2, 3, 5, 6]]
@test collect(celldofids(fespace)) == r
@test fespace.node_and_comp_to_dof ==  [-1, -2, -3, 1, 2, 3, 4, 5, 6]

trian = Triangulation(grid)
quad = CellQuadrature(trian,order=2)

bh = FEBasis(fespace)

uh = zero(fespace)

cellmat = integrate(inner(bh,bh),trian,quad)
cellvec = integrate(inner(bh,uh),trian,quad)

ufun(x) = x[1]

nfree = 6
ndiri = 3

test_fe_space(fespace, nfree, ndiri, cellmat, cellvec, ufun)

uh = interpolate(fespace,ufun)
u = CellField(trian,ufun)
e = u - uh

ufun1(x) = 1.0
ufun2(x) = 2.0
ufun5(x) = 5.0

dv = interpolate_diri_values(fespace,[ufun1,ufun2,ufun5])

@test dv == [1.0, 5.0, 2.0]

@test fespace.dirinodes == [1,2,3]

tol = 1.0e-8
@test sum(integrate(inner(e,e),trian,quad)) < tol

#writevtk(trian,"trian",nref=4,cellfields=["uh"=>uh])


# CLagrangianFESpace for vector values

T = VectorValue{2,Float64}

diritags = [1,2,5]
dirimasks = [(true,false),(false,true),(true,true)]

fespace = CLagrangianFESpace(
  T,grid,node_to_label,tag_to_labels,diritags,dirimasks)

ufun1(x) = VectorValue(1.0,-1.0)
ufun2(x) = VectorValue(2.0,-2.0)
ufun5(x) = VectorValue(5.0,-5.0)

dv = interpolate_diri_values(fespace,[ufun1,ufun2,ufun5])

@test dv == [1.0, 5.0, -5.0, -2.0]

ufun(x) = x

uh = interpolate(fespace,ufun)
u = CellField(trian,ufun)
e = u - uh

@test sum(integrate(inner(e,e),trian,quad)) < tol

bh = FEBasis(fespace)

uh = zero(fespace)

cellmat = integrate(inner(bh,bh),trian,quad)
cellvec = integrate(inner(bh,uh),trian,quad)

nfree = 14
ndiri = 4

test_fe_space(fespace, nfree, ndiri, cellmat, cellvec, ufun)

order = 1
fespace = CLagrangianFESpace(T,model,order,diritags,dirimasks)
test_fe_space(fespace, nfree, ndiri, cellmat, cellvec, ufun)

# Test the internal generation of high order meshes

model = CartesianDiscreteModel(partition=(1,2))
grid, node_to_label = grid_from_model_and_order(model,2)
points(grid)

r = Point{2,Float64}[
  (-1.0, -1.0), (1.0, -1.0), (-1.0, 0.0), (1.0, 0.0),
  (-1.0, 1.0), (1.0, 1.0), (0.0, -1.0), (0.0, 0.0),
  (-1.0, -0.5), (1.0, -0.5), (0.0, 1.0), (-1.0, 0.5),
  (1.0, 0.5), (0.0, -0.5), (0.0, 0.5)]

@test collect(points(grid)) ≈ r

r = [[1, 2, 3, 4, 7, 8, 9, 10, 14], [3, 4, 5, 6, 8, 11, 12, 13, 15]]

cells(grid)
@test collect(cells(grid)) == r

@test node_to_label == [1, 2, 7, 8, 3, 4, 5, 9, 7, 8, 6, 7, 8, 9, 9]

# Testing high order cases

model = CartesianDiscreteModel(partition=(1,2))

order = 2
T = VectorValue{2,Float64}
diritags = [1,2,5]
dirimasks = [(true,false),(false,true),(true,true)]

fespace = CLagrangianFESpace(T,model,order,diritags,dirimasks)

@test fespace.dirimasks === dirimasks

fespace = CLagrangianFESpace(T,model,order,diritags)

@test fespace.dirimasks == [(true, true), (true, true), (true, true)]

order = 2
T = Float64
diritags = ["physical_tag_1","physical_tag_2","physical_tag_5"]
dirimasks = [false,false,true]
fespace = CLagrangianFESpace(T,model,order,diritags,dirimasks)
@test fespace.dirimasks === dirimasks
@test fespace.diritags == [1,2,5]

fespace = CLagrangianFESpace(T,model,order,diritags)
@test fespace.dirimasks == [(true,), (true,), (true,)]

fespace = CLagrangianFESpace(T,model,order)
@test fespace.diritags == Int[]
@test fespace.dirimasks == Tuple{Bool}[]
@test all(fespace.node_and_comp_to_dof .> 0)
@test fespace.ndiridofs == 0

end # module
