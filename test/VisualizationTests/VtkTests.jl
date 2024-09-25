module VtkTests

using Test
using Gridap.Geometry: GridMock
using Gridap.Geometry: DiscreteModelMock
using Gridap.TensorValues
using Gridap.Arrays
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Visualization
using Gridap.Visualization: create_pvtk_file
using Gridap.CellData
using WriteVTK

d = mktempdir()
f = joinpath(d,"trian")

trian = GridMock()

node_ids = collect(1:num_nodes(trian))
cell_ids = collect(1:num_cells(trian))

mean(x) = sum(x)/length(x)

cell_center = lazy_map(mean, get_cell_coordinates(trian) )

for compress in [true,false]
  for append in [true,false]
    for ascii in [true,false]
      for vtkversion in [:default,:latest]
        write_vtk_file(
          trian,f,
          nodaldata=["nodeid"=>node_ids],
          celldata=["cellid"=>cell_ids,"centers"=>cell_center],
          compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
        )

        pvtk = Visualization.create_pvtk_file(
          trian,f; part=1, nparts=1,
          nodaldata=["nodeid"=>node_ids],
          celldata=["cellid"=>cell_ids,"centers"=>cell_center],
          compress=compress, append=append, ascii=ascii, vtkversion=vtkversion
        )
        vtk_save(pvtk)
      end
    end
  end
end

reffe = LagrangianRefFE(VectorValue{3,Float64},WEDGE,(3,3,4))
f = joinpath(d,"reffe")
writevtk(reffe,f)

reffe = LagrangianRefFE(VectorValue{2,Float64},QUAD,(2,0))
f = joinpath(d,"reffe")
writevtk(reffe,f)

f = joinpath(d,"poly")
writevtk(HEX,f)

domain = (0,1,0,1,0,1)
partition = (3,4,2)
grid = CartesianGrid(domain,partition)
model = UnstructuredDiscreteModel(grid)

f = joinpath(d,"model")
writevtk(model,f)

domain = (0,1,0,1,0,1)
partition = (3,4,2)
model = CartesianDiscreteModel(domain,partition)

f = joinpath(d,"model")
writevtk(model,f)

f = joinpath(d,"model")
model = DiscreteModelMock()
writevtk(model,f,labels=get_face_labeling(model))

f = joinpath(d,"trian")
trian = Triangulation(model)
writevtk(trian,f,order=2)

domain = (0,1,0,1)
partition = (2,4)
trian = Triangulation(CartesianDiscreteModel(domain,partition))

writevtk(trian,f,nsubcells=5,celldata=["rnd"=>rand(num_cells(trian))])

cf(x) = sin(4*x[1]*pi)*cos(5*x[2]*pi)

writevtk(trian,f,nsubcells=10, cellfields=[
  "cf"=>cf,
  "a"=>x->1,
  "v2"=>x->VectorValue(1,2),
  "v"=>x->VectorValue(1,2,3),
  "s"=>x->SymTensorValue(1.0,2.0,3.0),
  "c"=>x->SymFourthOrderTensorValue(1,2,3, 1,2,3, 1,2,3),
  "t"=>x->TensorValue(1,2,3,4),])

trian = GridMock()

p1 = Point{2,Float64}[(0.25,0.25),(0.75,0.75)]
p2 = Point{2,Float64}[(0.2,0.2),(0.4,0.4)]
q = CompressedArray([p1,p2], get_cell_type(trian))
q2x = get_cell_map(trian)
x = lazy_map(evaluate,q2x,q)

f = joinpath(d,"x")
writevtk(x,f,celldata=["cellid" => collect(1:num_cells(trian))], nodaldata = ["x" => x])


# Write VTK_LAGRANGE_* FE elements
writevtk(Grid(LagrangianRefFE(Float64,TRI,3)),joinpath(d,"tri_order3"))
writevtk(Grid(LagrangianRefFE(Float64,TRI,4)),joinpath(d,"tri_order4"))
writevtk(Grid(LagrangianRefFE(Float64,TRI,5)),joinpath(d,"tri_order5"))
writevtk(Grid(LagrangianRefFE(Float64,QUAD,3)),joinpath(d,"quad_order3"))
writevtk(Grid(LagrangianRefFE(Float64,QUAD,4)),joinpath(d,"quad_order4"))
writevtk(Grid(LagrangianRefFE(Float64,TET,3)),joinpath(d,"tet_order1"))
writevtk(Grid(LagrangianRefFE(Float64,HEX,3)),joinpath(d,"hex_order1"))

# Paraview collections
model = DiscreteModelMock()
trian = Triangulation(model)
f = joinpath(d,"collection")
paraview_collection(f) do pvd
    for i in 1:10
        pvd[Float64(i)] = createvtk(trian, f*"_$i", celldata=["rnd"=>rand(num_cells(trian))], cellfields=["cf" => cf])
        pvd[Float64(10+i)] = createvtk(x,f*"_$(10+i)",celldata=["cellid" => collect(1:num_cells(trian))], nodaldata = ["x" => x])
    end
    vtk_save(pvd)
end

# Paraview collections using createpvd
f = joinpath(d,"createpvd")
createpvd(f) do pvd
    for i in 1:10
        pvd[Float64(i)] = createvtk(trian, f*"_$i", celldata=["rnd"=>rand(num_cells(trian))], cellfields=["cf" => cf])
        pvd[Float64(10+i)] = createvtk(x,f*"_$(10+i)",celldata=["cellid" => collect(1:num_cells(trian))], nodaldata = ["x" => x])
    end
end
@test isfile(f*".pvd")

# Paraview collections using createpvd and parts
f = joinpath(d,"createpvd_parts")
parts = nothing
createpvd(parts,f) do pvd
    for i in 1:10
        pvd[Float64(i)] = createvtk(trian, f*"_$i", celldata=["rnd"=>rand(num_cells(trian))], cellfields=["cf" => cf])
        pvd[Float64(10+i)] = createvtk(x,f*"_$(10+i)",celldata=["cellid" => collect(1:num_cells(trian))], nodaldata = ["x" => x])
    end
end
@test isfile(f*".pvd")

f = joinpath(d,"x")
x = get_cell_points(CellQuadrature(trian,2))
writevtk(x,f;cellfields=["cf"=>cf])

## Visualize AppendedTriangulation
#
#domain = (0,1,0,1)
#partition = (10,10)
#grid1 = CartesianGrid(domain,partition)
#
#domain = (1,2,0,1)
#partition = (10,10)
#grid2 = simplexify(CartesianGrid(domain,partition))
#
#trian = lazy_append(grid1,grid2)
#
#f = joinpath(d,"trian")
#writevtk(trian,f)

f=joinpath(d,"pvdfile")
pvd = createpvd(f)
@test isa(pvd,WriteVTK.CollectionFile)
savepvd(pvd)
@test isfile(f*".pvd")

f=joinpath(d,"pvdfile_parts")
pvd = createpvd(nothing,f)
@test isa(pvd,WriteVTK.CollectionFile)
savepvd(pvd)
@test isfile(f*".pvd")

# Empty meshes
cells = (10,10)
domain = (0,1,0,1)
model = CartesianDiscreteModel(domain,cells)

# Domains and measures
Ω1 = Interior(model)
Ω2 = Interior(model,Int[])

f=joinpath(d,"empty")

write_vtk_file(Ω2,f,celldata=["u"=>rand(num_cells(Ω2))])
pvtk = create_pvtk_file(Ω1,f; part=1, nparts=2,celldata=["u"=>rand(num_cells(Ω1))])
vtk_save(pvtk)
pvtk = create_pvtk_file(Ω2,f; part=2, nparts=2,celldata=["u"=>rand(num_cells(Ω2))])
vtk_save(pvtk)


rm(d,recursive=true)

end # module
