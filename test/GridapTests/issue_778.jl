module Issue778

using Gridap

d = mktempdir()

domain = (0,1,0,1,0,1)
cells = (3,3,3)
model = CartesianDiscreteModel(domain,cells,isperiodic=(true,false,true))
model = simplexify(model)
writevtk(model,joinpath(d,"tmp_model_3"))

domain = (0,1,0,1)
cells = (3,3)
model = CartesianDiscreteModel(domain,cells,isperiodic=(true,false))
model = simplexify(model)
writevtk(model,joinpath(d,"tmp_model_2"))

rm(d,recursive=true)

end # module
