module DraftTests

using Gridap

partition = (0,1,0,1)
cells = (5,5)
model = CartesianDiscreteModel(partition,cells)

D = num_cell_dims(model)
Ω = Triangulation(ReferenceFE{D},model)
Γ = Triangulation(ReferenceFE{D-1},model)
∂K = CellBoundary(model)

v = CellField(x->x[1]+x[2],Ω)
q = CellField(x->x[1]-x[2]+1,Γ)

writevtk(Ω,"draft_Ω",cellfields=["v"=>v])
writevtk(Γ,"draft_Γ",cellfields=["q"=>q])
writevtk(∂K,"draft_∂K",offset=0.15,cellfields=["q"=>q,"v"=>v])

end # module
