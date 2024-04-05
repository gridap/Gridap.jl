using Gridap
using FillArrays
using Test

ptr  = [ 1, 5, 9, 13, 17 ]
data = [ 1,2,4,5, 2,3,5,6, 4,5,7,8, 5,6,8,9 ]
cell_vertex_lids = Gridap.Arrays.Table(data,ptr)
node_coordinates = Vector{Point{2,Float64}}(undef,11)

for j in 1:3
  for i in 1:3
    node_coordinates[3*(j-1)+i]=Point{2,Float64}((i-1)/2.0,(j-1)/2.0)
  end
end

polytope = QUAD
scalar_reffe = Gridap.ReferenceFEs.ReferenceFE(polytope,Gridap.ReferenceFEs.lagrangian,Float64,1)
cell_types = collect(Fill(1,length(cell_vertex_lids)))
cell_reffes = [scalar_reffe]
grid = Gridap.Geometry.UnstructuredGrid(node_coordinates,
                                        cell_vertex_lids,
                                        cell_reffes,
                                        cell_types,
                                        Gridap.Geometry.NonOriented())
model = Gridap.Geometry.UnstructuredDiscreteModel(grid)

Ω = Triangulation(model)
dΩ = Measure(Ω,2)

@test dΩ.quad.cell_point === dΩ.quad.cell_point # true
@test get_cell_points(dΩ) === get_cell_points(dΩ) # this was failing in the issue #905

# related application, taken from issue #905
order = 1
function project(order,Ω)
  degree = 2*order
  dΩ = Measure(Ω,degree)
  f = CellState(1.0,dΩ)
  reffe = ReferenceFE(lagrangian,Float64,order)
  V = FESpace(model,reffe,conformity=:L2)
  a(u,v) = ∫( u*v )dΩ
  l(v) = ∫( v*f )*dΩ
  op = AffineFEOperator(a,l,V,V)
  fh = solve(op)
  fh
end

# @test project(order,Ω) isa Any

@test try project(order,Ω)
  true
catch e
  rethrow(e)
  false
end
