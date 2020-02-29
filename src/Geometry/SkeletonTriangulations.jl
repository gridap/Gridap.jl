
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

"""
"""
function InterfaceTriangulation(model::DiscreteModel,cell_to_is_left::Vector{Bool})

  D = num_cell_dims(model)
  facet_grid = Grid(ReferenceFE{D-1},model)
  cell_grid = Grid(ReferenceFE{D},model)

  topo = get_grid_topology(model)
  facet_to_cells = Table(get_faces(topo,D-1,D))
  ifacet_to_facet, facet_to_lcell_left, facet_to_lcell_right = _find_interface_facets(
    cell_to_is_left, facet_to_cells)

  ifacet_trian = TriangulationPortion(facet_grid,ifacet_to_facet)

  left = GenericBoundaryTriangulation(ifacet_trian,cell_grid,topo,ifacet_to_facet,facet_to_lcell_left)
  right = GenericBoundaryTriangulation(ifacet_trian,cell_grid,topo,ifacet_to_facet,facet_to_lcell_right)

  SkeletonTriangulation(left,right)
end

function _find_interface_facets( cell_to_is_left, facet_to_cells::Table)

  nifacets = 0
  for facet in 1:length(facet_to_cells)
    a = facet_to_cells.ptrs[facet]
    b = facet_to_cells.ptrs[facet+1]
    if b-a == 2
      cell1 = facet_to_cells.data[a]
      cell2 = facet_to_cells.data[a+1]
      is_left_1 = cell_to_is_left[cell1]
      is_left_2 = cell_to_is_left[cell2]
      if is_left_1 != is_left_2
        nifacets += 1
      end
    end
  end

  T = eltype(eltype(facet_to_cells))
  ifacet_to_facet = zeros(T,nifacets)
  nfacets = length(facet_to_cells)
  facet_to_lcell_left = fill(Int8(1),nfacets)
  facet_to_lcell_right = fill(Int8(2),nfacets)

  nifacets = 0
  for facet in 1:length(facet_to_cells)
    a = facet_to_cells.ptrs[facet]
    b = facet_to_cells.ptrs[facet+1]
    if b-a == 2
      cell1 = facet_to_cells.data[a]
      cell2 = facet_to_cells.data[a+1]
      is_left_1 = cell_to_is_left[cell1]
      is_left_2 = cell_to_is_left[cell2]
      if is_left_1 != is_left_2
        nifacets += 1
        ifacet_to_facet[nifacets] = facet
        if is_left_1
          facet_to_lcell_left[facet] = 1
          facet_to_lcell_right[facet] = 2
        else
          facet_to_lcell_left[facet] = 2
          facet_to_lcell_right[facet] = 1
        end
      end
    end
  end

   ifacet_to_facet, facet_to_lcell_left, facet_to_lcell_right
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

function get_cell_id(trian::SkeletonTriangulation)
  left = get_cell_id(trian.left)
  right = get_cell_id(trian.right)
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

# Specific API

"""
"""
get_left_boundary(trian::SkeletonTriangulation) = trian.left

"""
"""
get_right_boundary(trian::SkeletonTriangulation) = trian.right

