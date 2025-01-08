
struct HybridTriangulation{Dc,Dp,A,B} <: Triangulation{Dc,Dp}
  cell_trian::A
  face_trian::B
  function HybridTriangulation(
    cell_trian :: BodyFittedTriangulation{Dc,Dp},
    face_trian :: BoundaryTriangulation{Df,Dp}
  ) where {Dc,Df,Dp}
    @assert Dc == Df+1 # TODO: Can we relax this?
    A, B = typeof(cell_trian), typeof(face_trian)
    new{Df,Dp,A,B}(cell_trian,face_trian)
  end
end

get_background_model(t::HybridTriangulation) = get_background_model(t.face_trian)
get_grid(t::HybridTriangulation) = get_grid(t.face_trian)
get_glue(t::HybridTriangulation,::Val{d}) where d = get_glue(t.face_trian,Val(d))
get_facet_normal(t::HybridTriangulation) = get_facet_normal(t.face_trian)

function HybridTriangulation(model::DiscreteModel{Dc},cell_to_bgcell) where Dc
  Df = Dc-1 # TODO: Expand to arbitrary face dimension

  # Cell triangulation
  cell_trian = Triangulation(model,cell_to_bgcell)

  # Get faces around each cell and their local index
  topo = get_grid_topology(model)
  bgcell_to_bgfaces = get_faces(topo,Dc,Df)
  bgface_to_bgcells = get_faces(topo,Df,Dc)

  nfaces = sum(bgcell -> length(view(bgcell_to_bgfaces,bgcell)), cell_to_bgcell)
  face_to_bgface = Vector{Int}(undef,nfaces)
  face_to_lcell = Vector{Int8}(undef,nfaces)

  face = 1
  for (cell,bgcell) in enumerate(cell_to_bgcell)
    bgfaces = view(bgcell_to_bgfaces,bgcell)
    for (lface,bgface) in enumerate(bgfaces)
      lcell = findfirst(x -> x == bgcell, view(bgface_to_bgcells,bgface))
      face_to_bgface[face] = bgface
      face_to_lcell[face] = lcell
      face += 1
    end
  end

  # Face triangulation around the selected cells
  cell_grid = get_grid(model)
  face_grid = view(Grid(ReferenceFE{Df},model),face_to_bgface)
  glue = OverlappingFaceToCellGlue(topo,cell_grid,face_grid,face_to_bgface,face_to_lcell)
  trian = BodyFittedTriangulation(model,face_grid,face_to_bgface)
  face_trian = BoundaryTriangulation(trian,glue)

  return HybridTriangulation(cell_trian,face_trian)
end
