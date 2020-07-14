"""
    struct PDiscRefFE{D} <: LagrangianRefFE{D}
      # Private fields
    end
"""
struct PDiscRefFE{D} <: LagrangianRefFE{D}
  reffe::LagrangianRefFE{D}
  polytope::Polytope{D}
end

function PDiscRefFE(::Type{T},p::Polytope,order::Integer) where T
  D = num_cell_dims(p)
  extrusion = tfill(TET_AXIS,Val{D}())
  simplex = ExtrusionPolytope(extrusion)
  reffe = LagrangianRefFE(T,simplex,order)
  PDiscRefFE{D}(reffe,p)
end

# LagrangianRefFE

get_face_nodes(reffe::PDiscRefFE) = get_face_own_nodes(reffe)

# Reffe

num_dofs(reffe::PDiscRefFE) = num_dofs(reffe.reffe)

get_polytope(reffe::PDiscRefFE) = reffe.polytope

get_prebasis(reffe::PDiscRefFE) = get_prebasis(reffe.reffe)

get_dof_basis(reffe::PDiscRefFE) = get_dof_basis(reffe.reffe)

get_default_conformity(reffe::PDiscRefFE) = L2Conformity()

get_face_dofs(reffe::PDiscRefFE) = get_face_own_dofs(reffe)

get_shapefuns(reffe::PDiscRefFE) = get_shapefuns(reffe.reffe)

