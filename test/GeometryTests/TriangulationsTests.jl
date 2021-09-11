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

@test model === get_discrete_model(Ω)
@test get_grid(Ω) === get_grid(model)
glue = get_glue(Ω,Val(2))
@test isa(glue.tface_to_mface,IdentityVector)
@test isa(glue.mface_to_tface,IdentityVector)
glue.mface_to_tface === glue.tface_to_mface
@test isa(glue.tface_to_mface_map,Fill)

Γ = Triangulation(ReferenceFE{1},model)
@test model === get_discrete_model(Γ)
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
glue = get_glue(Ω1,Val(2))
@test glue.tface_to_mface == findall(collect1d(cell_mask))
@test isa(glue.tface_to_mface_map,Fill)
@test model === get_discrete_model(Ω)
@test isa(glue.mface_to_tface,PosNegPartition)
@test glue.mface_to_tface[glue.tface_to_mface] == 1:length(glue.tface_to_mface)

cell_mcell = findall(collect1d(cell_mask))
Ω1 = Triangulation(model,cell_mcell)
glue = get_glue(Ω1,Val(2))
@test glue.tface_to_mface == cell_mcell
@test isa(glue.tface_to_mface_map,Fill)
@test model === get_discrete_model(Ω)
@test isa(glue.mface_to_tface,PosNegPartition)
@test glue.mface_to_tface[glue.tface_to_mface] == 1:length(glue.tface_to_mface)

Ω1 = view(Ω,cell_mcell)
glue = get_glue(Ω1,Val(2))
@test glue.tface_to_mface == cell_mcell
@test isa(glue.tface_to_mface_map,Fill)
@test model === get_discrete_model(Ω)
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
@test model === get_discrete_model(Ω)
@test isa(glue.mface_to_tface,PosNegPartition)
@test glue.mface_to_tface[glue.tface_to_mface] == 1:length(glue.tface_to_mface)

end # module
