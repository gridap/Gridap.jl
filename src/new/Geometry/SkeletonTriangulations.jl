
"""
    struct SkeletonTriangulation{Dc,Dp,B} <: Triangulation{Dc,Dp}
      left::B
      right::B
    end

The inner constructor enforces `B<:BoundaryTriangulation`
"""
struct SkeletonTriangulation{Dc,Dp,B} <: Triangulation{Dc,Dp}
  left::B
  right::B
  function SkeletonTriangulation(left::B,right::B) where B<:BoundaryTriangulation
    Dc = num_cell_dims(left)
    Dp = num_point_dims(left)
    new{Dc,Dp,B}(left,right)
  end
end

"""
    SkeletonTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
    SkeletonTriangulation(model::DiscreteModel)
"""
function SkeletonTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
  left_cell_around = 1
  left = GenericBoundaryTriangulation(model,face_to_mask,left_cell_around)
  right_cell_around = 2
  right = GenericBoundaryTriangulation(model,face_to_mask,right_cell_around)
  SkeletonTriangulation(left,right)
end

function SkeletonTriangulation(model::DiscreteModel)
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  face_to_mask = collect(Bool, .!get_isboundary_face(topo,D-1))
  SkeletonTriangulation(model,face_to_mask)
end

# Triangulation interface

function get_cell_coordinates(trian::SkeletonTriangulation)
  get_cell_coordinates(trian.left)
end

function get_reffes(trian::SkeletonTriangulation)
    get_reffes(trian.left)
end

function get_cell_type(trian::SkeletonTriangulation)
  get_cell_type(trian.left)
end

function reindex(a::AbstractArray,trian::SkeletonTriangulation)
  left = reindex(a,trian.left)
  right = reindex(a,trian.right)
  SkeletonPair(left,right)
end

function restrict(f::AbstractArray, trian::SkeletonTriangulation)
  left = restrict(f,trian.left)
  right = restrict(f,trian.right)
  SkeletonPair(left,right)
end

# Delegating into the left side

"""
    get_volume_triangulation(trian::SkeletonTriangulation)
"""
function get_volume_triangulation(trian::SkeletonTriangulation)
  get_volume_triangulation(trian.left)
end

"""
    get_normal_vector(trian::SkeletonTriangulation)
"""
function get_normal_vector(trian::SkeletonTriangulation)
  get_normal_vector(trian.left)
end
