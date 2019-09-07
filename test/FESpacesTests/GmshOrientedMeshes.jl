using Gridap

using GridapGmsh

model = GmshDiscreteModel("light-mesh/light-mesh.msh")

model.grids[2]
