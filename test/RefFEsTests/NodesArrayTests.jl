module NodesArrayTests

using Gridap, Test

p = Polytope(HEX_AXIS,HEX_AXIS)
order = (4,4)
vs = Gridap.Polytopes.generate_interior_nodes(p.nfaces[end], order)
vs_float = Gridap.Polytopes._equidistant_nodes_coordinates(vs,order)

p = Polytope(HEX_AXIS,TET_AXIS)
order = (3,3)
vs = Gridap.Polytopes.generate_interior_nodes(p.nfaces[end], order)
vs_float = Gridap.Polytopes._equidistant_nodes_coordinates(vs,order)

end #module
