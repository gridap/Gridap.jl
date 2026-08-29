module MTWRefFEsTest

using Test
using Gridap
using Gridap.ReferenceFEs
using Gridap.Geometry
using Gridap.Helpers
using Gridap.FESpaces
using Gridap.TensorValues
using Gridap.CellData
using Gridap.Arrays
using Gridap.Polynomials
using Gridap.Fields

using LinearAlgebra


partition = (1,1)
domain = (0.,1.,0.,1.)
model = simplexify(CartesianDiscreteModel(domain, partition))
#model = simplexify(CartesianDiscreteModel(domain, partition), positive=true)
#
mtw2 = MardalTaiWintherRefFE(Float64,TRI)
U = FESpace(model, mtw2)
#
Ω = Triangulation(model)
dΩ = Measure(Ω,8)
l2(u) = sum(∫(u⋅u)dΩ)
#
function φi_φih(i,U)
  #cφ(i) = collect(Float64(j==i) for j in 1:nfree )
  φic = collect(Float64(j==i) for j in 1:U.nfree )
  φi = FEFunction(U, φic)
  fφi = x -> φi(x)
  φih = interpolate(fφi, U)
  φi, φih
end
#
M = Matrix{Float64}(undef, (U.nfree,U.nfree));
for i in 1:U.nfree
  φi, φih = φi_φih(i,U)
  col = round.(φih.free_values; digits=3)
  setindex!(M, col, i, :)
  @show i, l2(φi-φih) ≤ 1.e-10
end
M

u(x) = VectorValue(0., 0.)
uh = interpolate(u, U);
l2(uh-u)

u(x) = VectorValue(1., 0.)
uh = interpolate(u, U);
l2(uh-u)

φ1, φ1h = φi_φih(1,U)
l2(φ1 -φ1h)
σ1 =


# works

function _polytope_model(p::Polytope{D}, v=nothing) where D
  model = UnstructuredDiscreteModel(UnstructuredGrid(ReferenceFE{D},p))
  if !isnothing(v)
    topo
  end
  label = get_face_labeling(model)
  for dface_to_entity in label.d_to_dface_to_entity[1:end-1] dface_to_entity .= 1 end
  add_tag!(label, "boundary", [1])
  model
end

model = _polytope_model(TRI)
#
setindex!(model.grid.node_coordinates, Point(1.,0), 1)
setindex!(model.grid.node_coordinates, Point(0, 1), 1)
setindex!(model.grid.node_coordinates, Point(1, 1), 1)
setindex!(model.grid_topology.vertex_coordinates, Point(1.,0), 1)
setindex!(model.grid_topology.vertex_coordinates, Point(0, 1), 1)
setindex!(model.grid_topology.vertex_coordinates, Point(1, 1), 1)
#
Ω_TRI = Triangulation(model)

U = FESpace(Ω_TRI, mtw2)
nfree = U.nfree
#
dΩ_TRI = Measure(Ω_TRI,6)
l2(u) = sum(∫(u⋅u )dΩ_TRI)
#
M = Matrix{Float64}(undef, (U.nfree,U.nfree));
for i in 1:nfree
  φi, φih = φi_φih(i,U)
  col = round.(φih.free_values; digits=3)
  setindex!(M, col, i, :)
  @show i, l2(φi-φih) ≤ 1.e-10
end
M


reffe = MardalTaiWintherRefFE(Float64,TET)

end
