module RestrictedDiscreteModelsTests

using Gridap.Arrays
using Gridap.ReferenceFEs
using Gridap.Geometry
using Test

domain = (0,2,0,2,0,2)
partition = (10,10,10)
oldmodel = CartesianDiscreteModel(domain,partition)

const R = 0.7

function is_in(coords)
  n = length(coords)
  x = (1/n)*sum(coords)
  d = x[1]^2 + x[2]^2 - R^2
  d < 0
end

oldgrid = get_grid(oldmodel)
oldcell_to_coods = get_cell_coordinates(oldgrid)
cell_to_mask = lazy_map(is_in,oldcell_to_coods)
cell_to_oldcell = findall(collect1d(cell_to_mask))

labels = get_face_labeling(oldmodel)

ne = num_entities(labels)
fluid_entity = ne+1
solid_entity = ne+2
cell_to_entity = get_cell_entity(labels)
cell_to_entity .= fluid_entity
cell_to_entity[cell_to_oldcell] .= solid_entity
add_tag!(labels,"fluid",[fluid_entity])
add_tag!(labels,"solid",[solid_entity])

model = RestrictedDiscreteModel(oldmodel,cell_to_oldcell)
test_discrete_model(model)

trian = get_triangulation(model)
@test isa(trian,RestrictedTriangulation)

trian = Triangulation(model)
@test isa(trian,RestrictedTriangulation)

model = DiscreteModel(oldmodel,cell_to_oldcell)
test_discrete_model(model)
@test isa(model,RestrictedDiscreteModel)

model = DiscreteModel(oldmodel,cell_to_mask)
test_discrete_model(model)
@test isa(model,RestrictedDiscreteModel)

model_fluid = DiscreteModel(oldmodel,tags="fluid")
test_discrete_model(model_fluid)
@test isa(model_fluid,RestrictedDiscreteModel)

model_solid = DiscreteModel(oldmodel,tags="solid")
test_discrete_model(model_solid)
@test isa(model_solid,RestrictedDiscreteModel)

itrian = InterfaceTriangulation(model_fluid,model_solid)

#using Gridap.Visualization
#writevtk(itrian,"itrian",nsubcells=3,cellfields=["normal"=>get_normal_vector(itrian)])

#
#writevtk(model,"model")
#writevtk(trian,"trian")

oldtrian = Triangulation(oldmodel)
trian_s = SkeletonTriangulation(model)
trian_b = BoundaryTriangulation(model)

@test get_background_triangulation(trian_s) === oldtrian
@test get_background_triangulation(trian_b) === oldtrian
@test get_background_triangulation(itrian) === oldtrian

#meas_K_b = cell_measure(trian_b,oldtrian)
#meas_K_sl = cell_measure(trian_s.plus,oldtrian)
#meas_K_sr = cell_measure(trian_s.minus,oldtrian)
#
#oldcell_to_cell = zeros(Int,num_cells(oldmodel))
#oldcell_to_cell[cell_to_oldcell] .= 1:length(cell_to_oldcell)
#
#@test all( oldcell_to_cell[findall(meas_K_b .!= 0)] .!= 0 )
#@test all( oldcell_to_cell[findall(meas_K_sl .!= 0)] .!= 0 )
#@test all( oldcell_to_cell[findall(meas_K_sr .!= 0)] .!= 0 )

#using Gridap.Visualization
#writevtk(trian_b,"trian_b")
#writevtk(trian_s,"trian_s")
#writevtk(trian,"trian")
#writevtk(oldtrian,"oldtrian",celldata=["meas_K_b"=>meas_K_b,"meas_K_sl"=>meas_K_sl,"meas_K_sr"=>meas_K_sr])

end # module
