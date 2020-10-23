module TriangulationsTests

using Test
using FillArrays
using LinearAlgebra: norm

using Gridap.Geometry
using Gridap.Arrays
using Gridap.TensorValues
using Gridap.Fields
using Gridap.ReferenceFEs
#using Gridap.CellData

import Gridap.Geometry: get_cell_type
import Gridap.Geometry: get_reffes
import Gridap.Geometry: get_cell_coordinates

struct MockTriangulation <: Triangulation{2,2} end

function get_cell_coordinates(::MockTriangulation)
  c1 = Point{2,Float64}[(0,0), (2,0), (0,2)]
  c2 = Point{2,Float64}[(2,0), (2,2), (1,1)]
  c3 = Point{2,Float64}[(2,2), (0,2), (1,1)]
  [c1,c2,c3]
end

function get_cell_type(::MockTriangulation)
  ncells = 3
  Fill(Int8(1),ncells)
end

function get_reffes(::MockTriangulation)
  [TRI3,]
end

trian = MockTriangulation()
test_triangulation(trian)

a = lazy_map(x->norm(x[1]-x[2]),get_cell_coordinates(trian))
test_array(a,collect(a))

cell_map = get_cell_map(trian)

ncells = num_cells(trian)
@test ncells == 3

qi = Point(0.5,0.5)
np = 4
qe = fill(qi,np)
q = Fill(qe,ncells)

xi1 = Point(1.0, 1.0)
xi2 = Point(1.5, 1.5)
xi3 = Point(0.5, 1.5)
x1 = fill(xi1,np)
x2 = fill(xi2,np)
x3 = fill(xi3,np)
x = [x1,x2,x3]

ji1 = TensorValue(2.0, 0.0, 0.0, 2.0)
ji2 = TensorValue(0.0, -1.0, 2.0, 1.0)
ji3 = TensorValue(-2.0, -1.0, 0.0, -1.0)
j1 = fill(ji1,np)
j2 = fill(ji2,np)
j3 = fill(ji3,np)
j = [j1,j2,j3]

test_array(cell_map,collect(cell_map))

lazy_map((ai,qi,vi,ji)->test_mapping(ai,(qi,),vi),cell_map,q,x)

lazy_map((ai,qi,vi,ji)->test_mapping(ai,(qi,),vi),cell_map,q,x)


lazy_map((ji,qi,vi)->test_mapping(ai,(qi,),vi),cell_map,q,x)

kk

test_field(cell_map,q,x,grad=j)

@test is_first_order(trian) == true

cf1 = CellField(3,trian)
cf2 = CellField(identity,trian)

x = get_physical_coordinate(trian)

@test get_cell_id(trian) == collect(1:num_cells(trian))
r = rand(num_cells(trian))
@test r === reindex(r,trian)

#using Gridap.Visualization
#writevtk(trian,"trian",cellfields=["cf1"=>cf1,"cf2"=>cf2,"x"=>x, "gradx"=>∇(x)])

end # module
