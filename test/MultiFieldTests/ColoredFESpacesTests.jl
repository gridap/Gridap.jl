using Gridap
using Gridap.Geometry, Gridap.FESpaces, Gridap.CellData, Gridap.ReferenceFEs, Gridap.Arrays, Gridap.Fields
using BlockArrays
using Gridap.MultiField

Dc = 2
model = CartesianDiscreteModel((0,1,0,1),(3,3))
topo = get_grid_topology(model)
labels = get_face_labeling(model)
add_tag_from_tags!(labels,"interface",["tag_1","tag_5","tag_2","tag_7","tag_3"])

Ω = Triangulation(model)

reffe = ReferenceFE(lagrangian,Float64,1)
V1 = FESpace(Ω,reffe;dirichlet_tags="tag_6")
V2 = FESpace(Ω,reffe;dirichlet_tags="tag_8")

V = MultiFieldFESpace([V1,V2])
W = MultiColorFESpace(V,[reffe,reffe],["interior","interface"],labels)

dΩ = Measure(Ω,4)

f = 1.0
a((u1,u2),(v1,v2)) = ∫(u1⋅v1 + u2⋅v2 - u1⋅v2)dΩ
l((v1,v2))   = ∫(f⋅v1)*dΩ

assem = SparseMatrixAssembler(W,W)
A = assemble_matrix(a,assem,W,W)
b = assemble_vector(l,assem,W)
