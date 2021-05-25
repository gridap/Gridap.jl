
"""
    struct SkeletonTriangulation{Dc,Dp,B} <: Triangulation{Dc,Dp}
      plus::B
      minus::B
    end
"""
struct SkeletonTriangulation{Dc,Dp,B} <: Triangulation{Dc,Dp}
  plus::B
  minus::B
  function SkeletonTriangulation(plus::B,minus::B) where B<:Triangulation
    Dc = num_cell_dims(plus)
    Dp = num_point_dims(plus)
    #@assert Dc + 1 == Dp
    new{Dc,Dp,B}(plus,minus)
  end
end

function Base.getproperty(x::SkeletonTriangulation, sym::Symbol)
  if sym == :⁺
    x.plus
  elseif sym == :⁻
    x.minus
  else
    getfield(x, sym)
  end
end

function Base.propertynames(x::SkeletonTriangulation, private::Bool=false)
  (fieldnames(typeof(x))...,:⁺,:⁻)
end

#have_compatible_domains(a::SkeletonTriangulation,b::Triangulation) = a.plus===b || a.minus===b
#have_compatible_domains(a::Triangulation,b::SkeletonTriangulation) = have_compatible_domains(b,a)
#have_compatible_domains(a::SkeletonTriangulation,b::SkeletonTriangulation) = a===b

"""
    SkeletonTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
    SkeletonTriangulation(model::DiscreteModel)
"""
function SkeletonTriangulation(model::DiscreteModel,face_to_mask::AbstractVector{Bool})
  left_cell_around = 1
  plus = BoundaryTriangulation(model,face_to_mask,left_cell_around)
  right_cell_around = 2
  minus = BoundaryTriangulation(model,face_to_mask,right_cell_around)
  SkeletonTriangulation(plus,minus)
end

function SkeletonTriangulation(model::DiscreteModel)
  topo = get_grid_topology(model)
  D = num_cell_dims(model)
  face_to_mask = collect(Bool, .!get_isboundary_face(topo,D-1))
  SkeletonTriangulation(model,face_to_mask)
end

function SkeletonTriangulation(
  model::DiscreteModel,
  reffe::ReferenceFE,
  face_own_dofs::Vector{Vector{Int}})

  ncells = num_cells(model)
  cell_to_reffe = Fill(1,ncells)
  SkeletonTriangulation(model,[reffe],cell_to_reffe,[face_own_dofs])
end

function SkeletonTriangulation(
  model::DiscreteModel,
  reffes::Vector{<:ReferenceFE},
  cell_to_reffe::AbstractVector,
  reffe_to_face_own_dofs::Vector{Vector{Vector{Int}}})

  reffe_to_polytope = map(get_polytope,reffes)
  reffe_to_faces = map(get_faces,reffe_to_polytope)
  reffe_to_range = map(p->get_dimrange(p,num_dims(p)-1),reffe_to_polytope)

  reffe_to_lfacet_to_mask = _compute_reffe_to_facet_to_mask(
    reffe_to_faces,reffe_to_range,reffe_to_face_own_dofs)

  topo = get_grid_topology(model)

  D = num_cell_dims(model)
  nfacets = num_facets(model)
  facet_to_isboundary = get_isboundary_face(topo,D-1)
  cell_to_facets = get_faces(topo,D,D-1)
  facet_to_mask = _compute_disc_facet_mask(
    facet_to_isboundary,cell_to_facets,cell_to_reffe,reffe_to_lfacet_to_mask)

  SkeletonTriangulation(model,facet_to_mask)

end

function _compute_disc_facet_mask(
  facet_to_isboundary, cell_to_facets::Table, cell_to_reffe, reffe_to_lfacet_to_mask)
  nfacets = length(facet_to_isboundary)
  facet_to_mask = Vector{Bool}(undef,nfacets)
  ncells = length(cell_to_reffe)
  for cell in 1:ncells
    pini = cell_to_facets.ptrs[cell]
    pend = cell_to_facets.ptrs[cell+1]-1
    reffe = cell_to_reffe[cell]
    lfacet_to_mask = reffe_to_lfacet_to_mask[reffe]
    for (lfacet,p) in enumerate(pini:pend)
      facet = cell_to_facets.data[p]
      mask = lfacet_to_mask[lfacet]
      isinterior = ! facet_to_isboundary[facet]
      facet_to_mask[facet] = mask && isinterior
    end
  end
  facet_to_mask
end

function _compute_reffe_to_facet_to_mask(
  reffe_to_faces,reffe_to_range,reffe_to_face_own_dofs)

  reffe_to_facet_to_mask = Vector{Bool}[]
  nreffes = length(reffe_to_faces)
  for reffe in 1:nreffes
    range = reffe_to_range[reffe]
    facet_to_faces = reffe_to_faces[reffe][range]
    facet_to_own_dofs = reffe_to_face_own_dofs[reffe][range]
    face_to_own_dofs = reffe_to_face_own_dofs[reffe]
    nfacets = length(facet_to_faces)
    facet_to_mask = fill(false,nfacets)
    for facet in 1:nfacets
      n = length(facet_to_own_dofs[facet])
      faces = facet_to_faces[facet]
      for face in faces
        n += length(face_to_own_dofs[face])
      end
      facet_to_mask[facet] = n == 0
    end
    push!(reffe_to_facet_to_mask,facet_to_mask)
  end

  reffe_to_facet_to_mask
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

  ifacet_trian = RestrictedTriangulation(facet_grid,ifacet_to_facet)

  glue_left = FaceToCellGlue(topo,cell_grid,ifacet_trian,ifacet_to_facet,facet_to_lcell_left)
  glue_right = FaceToCellGlue(topo,cell_grid,ifacet_trian,ifacet_to_facet,facet_to_lcell_right)

  plus = BoundaryTriangulation(ifacet_trian,cell_grid,glue_left)
  minus = BoundaryTriangulation(ifacet_trian,cell_grid,glue_right)

  SkeletonTriangulation(plus,minus)
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
  get_cell_coordinates(trian.plus)
end

function get_node_coordinates(trian::SkeletonTriangulation)
  get_node_coordinates(trian.plus)
end

function get_cell_node_ids(trian::SkeletonTriangulation)
  get_cell_node_ids(trian.plus)
end

function get_reffes(trian::SkeletonTriangulation)
    get_reffes(trian.plus)
end

function get_cell_type(trian::SkeletonTriangulation)
  get_cell_type(trian.plus)
end

function get_cell_map(trian::SkeletonTriangulation)
  get_cell_map(trian.plus)
end

function get_facet_normal(trian::SkeletonTriangulation)
  plus = get_facet_normal(trian.plus)
  #minus = get_facet_normal(trian.minus)
  minus = lazy_map(Broadcasting(Operation(-)),plus)
  SkeletonPair(plus,minus)
end

TriangulationStyle(::Type{<:SkeletonTriangulation}) = SubTriangulation()

function get_background_triangulation(trian::SkeletonTriangulation)
  get_background_triangulation(trian.plus)
end

function get_cell_to_bgcell(trian::SkeletonTriangulation)
  plus = get_cell_to_bgcell(trian.plus)
  minus = get_cell_to_bgcell(trian.minus)
  SkeletonPair(plus,minus)
end

function get_cell_ref_map(trian::SkeletonTriangulation)
  plus = get_cell_ref_map(trian.plus)
  minus = get_cell_ref_map(trian.minus)
  SkeletonPair(plus,minus)
end

function RestrictedTriangulation(
  oldtrian::SkeletonTriangulation,cell_to_oldcell::AbstractVector{<:Integer})
  plus = RestrictedTriangulation(oldtrian.plus,cell_to_oldcell)
  minus = RestrictedTriangulation(oldtrian.minus,cell_to_oldcell)
  SkeletonTriangulation(plus,minus)
end
