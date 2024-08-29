module TriangulationsTests

using Gridap.Geometry
using Gridap.Arrays
using Gridap.ReferenceFEs
using Test
using FillArrays

domain = (-1,1,-1,1)
cells = (10,10)
model = CartesianDiscreteModel(domain,cells)
Ω = Triangulation(model)
test_triangulation(Ω)
Ω_ = Triangulation(Ω)
@test Ω_ === Ω

@test model === get_background_model(Ω)
@test model === get_active_model(Ω)
@test get_grid(Ω) === get_grid(model)
glue = get_glue(Ω,Val(2))
@test isa(glue.tface_to_mface,IdentityVector)
@test isa(glue.mface_to_tface,IdentityVector)
glue.mface_to_tface === glue.tface_to_mface
@test isa(glue.tface_to_mface_map,Fill)

glue_0 = get_glue(Ω,Val(0))
glue_1 = get_glue(Ω,Val(1))
glue_2 = get_glue(Ω,Val(2))

Θ = GenericTriangulation(
  get_grid(Ω),
  get_background_model(Ω),
  (glue_0,glue_1,glue_2))
test_triangulation(Θ)
@test get_glue(Θ,Val(0)) === glue_0
@test get_glue(Θ,Val(1)) === glue_1
@test get_glue(Θ,Val(2)) === glue_2

Θ = GenericTriangulation(
  get_grid(Ω),
  get_background_model(Ω))
test_triangulation(Θ)

Θ = GenericTriangulation(get_grid(Ω))
@test isa(Θ,Triangulation)

Γ = Triangulation(ReferenceFE{1},model)
@test model === get_background_model(Γ)
glue = get_glue(Γ,Val(1))
@test isa(glue.tface_to_mface,IdentityVector)
@test isa(glue.mface_to_tface,IdentityVector)
glue.mface_to_tface === glue.tface_to_mface
@test isa(glue.tface_to_mface_map,Fill)

cell_xs = get_cell_coordinates(Ω)
cell_mask = lazy_map(cell_xs) do xs
  R = 0.7
  n = length(xs)
  x = (1/n)*sum(xs)
  d = x[1]^2 + x[2]^2 - R^2
  d < 0
end

Ω1 = Triangulation(model,cell_mask)
@test model === get_background_model(Ω1)
@test model !== get_active_model(Ω1)
glue = get_glue(Ω1,Val(2))
@test glue.tface_to_mface == findall(collect1d(cell_mask))
@test isa(glue.tface_to_mface_map,Fill)
@test model === get_background_model(Ω)
@test isa(glue.mface_to_tface,PosNegPartition)
@test glue.mface_to_tface[glue.tface_to_mface] == 1:length(glue.tface_to_mface)

tface_to_val = [ rand(3,4) for i in 1:num_cells(Ω1) ]
mface_to_val = extend(tface_to_val,glue.mface_to_tface)
@test mface_to_val[glue.tface_to_mface] == tface_to_val

cell_mcell = findall(collect1d(cell_mask))
Ω1 = Triangulation(model,cell_mcell)
glue = get_glue(Ω1,Val(2))
@test glue.tface_to_mface == cell_mcell
@test isa(glue.tface_to_mface_map,Fill)
@test model === get_background_model(Ω)
@test isa(glue.mface_to_tface,PosNegPartition)
@test glue.mface_to_tface[glue.tface_to_mface] == 1:length(glue.tface_to_mface)

Ω1 = view(Ω,cell_mcell)
glue = get_glue(Ω1,Val(2))
@test glue.tface_to_mface == cell_mcell
@test isa(glue.tface_to_mface_map,Fill)
@test model === get_background_model(Ω)
@test isa(glue.mface_to_tface,PosNegPartition)
@test glue.mface_to_tface[glue.tface_to_mface] == 1:length(glue.tface_to_mface)

labels = get_face_labeling(model)
entity = num_entities(labels)+1
labels.d_to_dface_to_entity[end][cell_mcell] .= entity
add_tag!(labels,"Ω1",[entity])
Ω1 = Triangulation(model,tags="Ω1")
glue = get_glue(Ω1,Val(2))
@test glue.tface_to_mface == cell_mcell
@test isa(glue.tface_to_mface_map,Fill)
@test model === get_background_model(Ω)
@test isa(glue.mface_to_tface,PosNegPartition)
@test glue.mface_to_tface[glue.tface_to_mface] == 1:length(glue.tface_to_mface)

glue = get_glue(Ω,Val(2))
tface_to_data = rand(num_cells(Ω))
mface_to_data = extend(tface_to_data,glue.mface_to_tface)
@test tface_to_data === mface_to_data

glue = get_glue(Ω1,Val(2))
tface_to_data = rand(num_cells(Ω1))
mface_to_data = extend(tface_to_data,glue.mface_to_tface)
@test mface_to_data[glue.tface_to_mface] == tface_to_data
tface_to_data = get_cell_shapefuns(Ω1)
mface_to_data = extend(tface_to_data,glue.mface_to_tface)
@test mface_to_data[glue.tface_to_mface] == tface_to_data

@test is_change_possible(Ω,Ω1)
@test is_change_possible(Ω1,Ω)

Ω2 = best_target(Ω1,Ω)
glue = get_glue(Ω2,Val(2))
@test isa(glue.tface_to_mface,IdentityVector)
@test isa(glue.mface_to_tface,IdentityVector)

# Using a non-injective tface_to_mface map
Ω3 = Triangulation(model, [1,2,1,2,3,3])
glue = get_glue(Ω3,Val(2))
@test isa(glue.tface_to_mface_map,Fill)
@test isa(glue.mface_to_tface,Nothing)

end # module
