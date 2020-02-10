module ExtrusionPolytopesTests

include("../../src/ReferenceFEs/ReferenceFEs.jl")

using Gridap.Fields
using .ReferenceFEs

using .ReferenceFEs: _polytopenfaces
using .ReferenceFEs: NFace
using .ReferenceFEs: DFace
using .ReferenceFEs: _polytopemesh
using .ReferenceFEs: _dimfrom_fs_dimto_fs
using .ReferenceFEs: _nface_to_nfacedim
using .ReferenceFEs: _nfaces_vertices
using .ReferenceFEs: _face_normals
using .ReferenceFEs: _edge_tangents
using .ReferenceFEs: _admissible_permutations
using .ReferenceFEs: _admissible_permutations_n_cube

extrusion = Point{0,Int}()
p = DFace(extrusion)
perms = _admissible_permutations(p)
display(perms)

extrusion = Point(HEX_AXIS)
p = DFace(extrusion)
perms = _admissible_permutations(p)
display(perms)

extrusion = Point(TET_AXIS)
p = DFace(extrusion)
perms = _admissible_permutations(p)
display(perms)

extrusion = Point(TET_AXIS,TET_AXIS)
p = DFace(extrusion)
perms = _admissible_permutations(p)
display(perms)

extrusion = Point(HEX_AXIS,HEX_AXIS)
p = DFace(extrusion)
perms = _admissible_permutations(p)
display(perms)

extrusion = Point(TET_AXIS, TET_AXIS, TET_AXIS)
p = DFace(extrusion)
perms = _admissible_permutations(p)
display(perms)

extrusion = Point(HEX_AXIS, HEX_AXIS, HEX_AXIS)
p = DFace(extrusion)
perms = _admissible_permutations(p)
display(perms)


f = DFace{2}(p,21)

f = DFace{3}(p,27)


end # module
