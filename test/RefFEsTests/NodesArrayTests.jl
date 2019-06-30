module NodesArrayTests

using Gridap, Test

p = Polytope(HEX_AXIS,HEX_AXIS)
vcs = vertices_coordinates(p)
@show vcs
p = Polytope(HEX_AXIS,TET_AXIS)
vcs = vertices_coordinates(p)
@show vcs


function vertices_coordinates(p::Polytope)
  vs = Gridap.Polytopes._dimfrom_fs_dimto_fs(p, dim, 0)[1]
  vcs = Point{dim,Float64}[]
  for i in 1:length(vs)
    cs = convert(Vector{Float64}, [p.nfaces[vs[i]].anchor...])
    push!(vcs,cs)
  end
  return vcs
end


end #module
