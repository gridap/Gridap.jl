
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
  function SkeletonTriangulation(left::B,right::B) where B<:Triangulation
    Dc = num_cell_dims(left)
    Dp = num_point_dims(left)
    @assert Dc + 1 == Dp
    new{Dc,Dp,B}(left,right)
  end
end

"""
    SkeletonTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
    SkeletonTriangulation(model::DiscreteModel)
"""
function SkeletonTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
  left_cell_around = 1
  left = BoundaryTriangulation(model,face_to_mask,left_cell_around)
  right_cell_around = 2
  right = BoundaryTriangulation(model,face_to_mask,right_cell_around)
  SkeletonTriangulation(left,right)
end

function SkeletonTriangulation(model::DiscreteModel)
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  face_to_mask = collect(Bool, .!get_isboundary_face(topo,D-1))
  SkeletonTriangulation(model,face_to_mask)
end

const IN = -1
const OUT = 1

"""
"""
function InterfaceTriangulation(model::DiscreteModel,cell_to_is_in::Vector{Bool})
  cell_to_inout = fill(Int8(OUT),length(cell_to_is_in))
  cell_to_inout[cell_to_is_in] .= IN
  InterfaceTriangulation(model,cell_to_inout)
end

function InterfaceTriangulation(model::DiscreteModel,cells_in,cells_out)
  cell_to_inout = fill(Int8(OUT),num_cells(model))
  cell_to_inout[cells_in] .= IN
  cell_to_inout[cells_out] .= OUT
  InterfaceTriangulation(model,cell_to_inout)
end

function InterfaceTriangulation(model::DiscreteModel,cell_to_inout::AbstractVector{<:Integer})

  D = num_cell_dims(model)
  facet_grid = Grid(ReferenceFE{D-1},model)
  cell_grid = Grid(ReferenceFE{D},model)

  topo = get_grid_topology(model)
  facet_to_cells = Table(get_faces(topo,D-1,D))
  ifacet_to_facet, facet_to_lcell_left, facet_to_lcell_right = _find_interface_facets(
    cell_to_inout, facet_to_cells)

  ifacet_trian = TriangulationPortion(facet_grid,ifacet_to_facet)

  left = GenericBoundaryTriangulation(ifacet_trian,cell_grid,topo,ifacet_to_facet,facet_to_lcell_left)
  right = GenericBoundaryTriangulation(ifacet_trian,cell_grid,topo,ifacet_to_facet,facet_to_lcell_right)

  SkeletonTriangulation(left,right)
end

function _find_interface_facets( cell_to_inout, facet_to_cells::Table)

  nifacets = 0
  for facet in 1:length(facet_to_cells)
    a = facet_to_cells.ptrs[facet]
    b = facet_to_cells.ptrs[facet+1]
    if b-a == 2
      cell1 = facet_to_cells.data[a]
      cell2 = facet_to_cells.data[a+1]
      inout_1 = cell_to_inout[cell1]
      inout_2 = cell_to_inout[cell2]
      if (inout_1== IN && inout_2==OUT) || (inout_1== OUT && inout_2==IN)
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
      inout_1 = cell_to_inout[cell1]
      inout_2 = cell_to_inout[cell2]
      if (inout_1== IN && inout_2==OUT) || (inout_1== OUT && inout_2==IN)
        nifacets += 1
        ifacet_to_facet[nifacets] = facet
        if inout_1 == IN
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

function TriangulationPortion(oldtrian::SkeletonTriangulation,cell_to_oldcell::Vector{Int})
  left = TriangulationPortion(oldtrian.left,cell_to_oldcell)
  right = TriangulationPortion(oldtrian.right,cell_to_oldcell)
  SkeletonTriangulation(left,right)
end

# Specific API

"""
"""
get_left_boundary(trian::SkeletonTriangulation) = trian.left

"""
"""
get_right_boundary(trian::SkeletonTriangulation) = trian.right

