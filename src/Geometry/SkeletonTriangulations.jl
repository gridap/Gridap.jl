"""
    struct SkeletonPair{L,R} <: GridapType
      plus::L
      minus::R
    end
"""
struct SkeletonPair{L,R} <: GridapType
  plus::L
  minus::R
end

function Base.getproperty(x::SkeletonPair, sym::Symbol)
  if sym == :⁺
    x.plus
  elseif sym == :⁻
    x.minus
  else
    getfield(x, sym)
  end
end

function Base.propertynames(x::SkeletonPair, private::Bool=false)
  (fieldnames(typeof(x))...,:⁺,:⁻)
end

"""
    struct SkeletonTriangulation{Dc,Dp,B} <: Triangulation{Dc,Dp}
      plus::B
      minus::B
    end
"""
struct SkeletonTriangulation{Dc,Dp,B,C} <: Triangulation{Dc,Dp}
  plus::B
  minus::C
  function SkeletonTriangulation(plus::B,minus::C) where {B<:Triangulation,C<:Triangulation}
    DcP = num_cell_dims(plus)
    DpP = num_point_dims(plus)
    DcM = num_cell_dims(minus)
    DpP = num_point_dims(minus)
    @assert DcP == DcM
    @assert DpP == DpP
    #@assert Dc + 1 == Dp
    new{DcP,DpP,B,C}(plus,minus)
  end
end

"""
    Skeleton(args...; kwargs...)

Alias for [`SkeletonTriangulation`](@ref)(args..., kwargs...).
"""
function Skeleton(args...;kwargs...)
  SkeletonTriangulation(args...;kwargs...)
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

# Triangulation interface

get_background_model(t::SkeletonTriangulation) = get_background_model(t.plus)
get_grid(t::SkeletonTriangulation) = get_grid(t.plus)
get_glue(t::SkeletonTriangulation{D},::Val{D}) where D = get_glue(t.plus,Val(D))

function get_glue(trian::SkeletonTriangulation,::Val{Dp}) where Dp
  model = get_background_model(trian)
  Dm = num_cell_dims(model)
  get_glue(trian,Val(Dp),Val(Dm))
end

function get_glue(trian::SkeletonTriangulation,::Val{Dp},::Val{Dm}) where {Dp,Dm}
  nothing
end

function get_glue(trian::SkeletonTriangulation,::Val{D},::Val{D}) where D
  plus = get_glue(trian.plus,Val(D))
  minus = get_glue(trian.minus,Val(D))
  SkeletonPair(plus,minus)
end

function is_change_possible(sglue::FaceToFaceGlue,tglue::SkeletonPair)
  is_change_possible(sglue,tglue.plus) && is_change_possible(sglue,tglue.minus)
end

function get_facet_normal(trian::SkeletonTriangulation)
  plus = get_facet_normal(trian.plus)
  minus = get_facet_normal(trian.minus)
  SkeletonPair(plus,minus)
end

function get_edge_tangent(trian::SkeletonTriangulation)
  plus = get_edge_tangent(trian.plus)
  minus = get_edge_tangent(trian.minus)
  SkeletonPair(plus,minus)
end

# Related with CompositeTriangulation
function _compose_glues(rglue::FaceToFaceGlue,dglue::SkeletonPair)
  plus = _compose_glues(rglue,dglue.plus)
  minus = _compose_glues(rglue,dglue.minus)
  SkeletonPair(plus,minus)
end

function Base.view(glue::SkeletonPair{<:FaceToFaceGlue},ids::AbstractArray)
  SkeletonPair(view(glue.plus,ids),view(glue.minus,ids))
end

function restrict(a::SkeletonPair,b::AbstractArray)
  SkeletonPair(restrict(a.plus,b),restrict(a.minus,b))
end

# Constructors

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

function SkeletonTriangulation(rtrian::Triangulation,args...;kwargs...)
  rmodel = get_active_model(rtrian)
  dtrian = SkeletonTriangulation(rmodel,args...;kwargs...)
  CompositeTriangulation(rtrian,dtrian)
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

"""
    Interface(args...; kwargs...)

Alias for [`InterfaceTriangulation`](@ref)(args...; kwargs...).
"""
function Interface(args...;kwargs...)
  InterfaceTriangulation(args...;kwargs...)
end

function InterfaceTriangulation(model::DiscreteModel,cells_in,cells_out)
  cell_to_inout = fill(Int8(0),num_cells(model))
  cell_to_inout[cells_in] .= IN
  cell_to_inout[cells_out] .= OUT
  InterfaceTriangulation(model,cell_to_inout)
end

function InterfaceTriangulation(trian_in::Triangulation,trian_out::Triangulation)
  d = num_cell_dims(trian_in)
  @assert d == num_cell_dims(trian_out)
  model = get_background_model(trian_in)
  D = num_cell_dims(model)
  if d == D
    _interface_between_volumes(trian_in,trian_out)
  elseif d+1 == D
    _interface_between_surfaces(trian_in,trian_out)
  else
    error("Not implemented")
  end
end

function _interface_between_volumes(trian_in,trian_out)
  D = num_cell_dims(trian_in)
  @assert D == num_cell_dims(trian_out)
  glue_in = get_glue(trian_in,Val(D))
  glue_out = get_glue(trian_out,Val(D))
  cells_in = glue_in.tface_to_mface
  cells_out = glue_out.tface_to_mface
  model = get_background_model(trian_in)
  @check model === get_background_model(trian_out)
  @notimplementedif D != num_cell_dims(model) "Not implemented, but it should be easy to implement."
  InterfaceTriangulation(model,cells_in,cells_out)
end

function _interface_between_surfaces(trian_in,trian_out)
  d = num_cell_dims(trian_in)
  @assert d == num_cell_dims(trian_out)
  glue_in = get_glue(trian_in,Val(d))
  glue_out = get_glue(trian_out,Val(d))
  faces_in = glue_in.tface_to_mface
  faces_out = glue_out.tface_to_mface
  faces = vcat(faces_in,faces_out)
  model = get_background_model(trian_in)
  D = num_cell_dims(model)
  @assert D == d+1 "Not implemented"
  rtrian = Boundary(model,faces) # TODO we take the first local cell around
  rmodel = get_active_model(rtrian)
  rfaces_in = collect(Int32,1:length(faces_in))
  rfaces_out = collect(Int32,length(faces_in) .+ (1:length(faces_out)))
  rtrian_in = Triangulation(rmodel,rfaces_in)
  rtrian_out = Triangulation(rmodel,rfaces_out)
  dtrian = InterfaceTriangulation(rtrian_in,rtrian_out)
  CompositeTriangulation(rtrian,dtrian)
end

function InterfaceTriangulation(model::DiscreteModel,cell_to_inout::AbstractVector{<:Integer})

  D = num_cell_dims(model)
  facet_grid = Grid(ReferenceFE{D-1},model)
  cell_grid = Grid(ReferenceFE{D},model)

  topo = get_grid_topology(model)
  facet_to_cells = Table(get_faces(topo,D-1,D))
  ifacet_to_facet, facet_to_lcell_left, facet_to_lcell_right = _find_interface_facets(
    cell_to_inout, facet_to_cells)

  ifacet_grid = view(facet_grid,ifacet_to_facet)

  glue_left = FaceToCellGlue(topo,cell_grid,ifacet_grid,ifacet_to_facet,facet_to_lcell_left)
  glue_right = FaceToCellGlue(topo,cell_grid,ifacet_grid,ifacet_to_facet,facet_to_lcell_right)

  trian = BodyFittedTriangulation(model,ifacet_grid,ifacet_to_facet)
  plus = BoundaryTriangulation(trian,glue_left)
  minus = BoundaryTriangulation(trian,glue_right)

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

