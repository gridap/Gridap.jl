
module HDivRefFEs

##

using Gridap, Test


p = Polytope(1,1)
d = dim(p)
facets = nfaces_dim(p,d-1)
verts = vertices_coordinates(p)
facet_vs = nface_connections(p,d-1,0)
##

facet_ps = nface_ref_polytopes(p)[facets]
@assert(all(extrusion.(facet_ps) .== extrusion(facet_ps[1])),
        "All facet must be of the same type")




end # module HDivRefFEs
