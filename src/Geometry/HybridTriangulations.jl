
struct CellToFaceGlue{A,B,C} <: GridapType
  tcell_to_mcell::A
  mcell_to_mfaces::B
  face_glue::C
end

function CellToFaceGlue(
  model::DiscreteModel{Dc},
  tcell_to_mcell,
  Df = Dc-1
) where Dc
  topo = get_grid_topology(model)
  mcell_to_mfaces = get_faces(topo,Dc,Df)
  mface_to_mcells = get_faces(topo,Df,Dc)

  tcell_to_mfaces = Table(mcell_to_mfaces[tcell_to_mcell])
  tface_to_mface = tcell_to_mfaces.data
  nfaces = length(tface_to_mface)

  tface_to_tcell = flatten_partition(Table(Base.OneTo(nfaces),tcell_to_mfaces.ptrs))
  tface_to_mcell = view(tcell_to_mcell,tface_to_tcell)
  tface_to_lcell = Arrays.find_local_nbor_index(tface_to_mcell, tface_to_mface, mface_to_mcells)

  cell_grid = get_grid(model)
  face_grid = view(Grid(ReferenceFE{Df},model),tface_to_mface)
  face_glue = OverlappingFaceToCellGlue(topo,cell_grid,face_grid,tface_to_mface,tface_to_lcell)
  return CellToFaceGlue(tcell_to_mcell,mcell_to_mfaces,face_glue)
end

function OverlappingFaceToCellGlue(
  topo::GridTopology,
  cell_grid::Grid,
  face_grid::Grid,
  face_to_bgface::AbstractVector,
  face_to_lcell::AbstractVector
)
  Dc = num_cell_dims(cell_grid)
  Df = num_cell_dims(face_grid)
  bgface_to_cells = get_faces(topo,Df,Dc)
  cell_to_bgfaces = get_faces(topo,Dc,Df)
  cell_to_lface_to_pindex = Table(get_cell_permutations(topo,Df))

  face_to_cells = lazy_map(Reindex(bgface_to_cells), face_to_bgface)
  face_to_cell  = collect(Int32,lazy_map(getindex,face_to_cells,face_to_lcell))
  face_to_lface = find_local_index(face_to_bgface,face_to_cell,cell_to_bgfaces)

  f = (p) -> fill(Int8(UNSET), num_faces(p,Df))
  ctype_to_lface_to_ftype = map(f, get_polytopes(cell_grid))
  face_to_ftype = get_cell_type(face_grid)
  cell_to_ctype = get_cell_type(cell_grid)

  _fill_ctype_to_lface_to_ftype!(
    ctype_to_lface_to_ftype,
    face_to_cell,
    face_to_lface,
    face_to_ftype,
    cell_to_ctype
  )

  FaceToCellGlue(
    face_to_bgface,
    face_to_cell,
    face_to_lface,
    face_to_lcell,
    face_to_ftype,
    cell_to_ctype,
    cell_to_lface_to_pindex,
    ctype_to_lface_to_ftype
  )
end


struct HybridTriangulation{Dc,Dp,A,B,C} <: Triangulation{Dc,Dp}
  cell_trian::A
  face_trian::B
  glue::C
  function HybridTriangulation(
    cell_trian :: BodyFittedTriangulation{Dc,Dp},
    face_trian :: BodyFittedTriangulation{Df,Dp},
    glue       :: CellToFaceGlue
  ) where {Dc,Df,Dp}
    @check get_background_model(cell_trian) === get_background_model(face_trian)
    A, B, C = typeof(cell_trian), typeof(face_trian), typeof(glue)
    new{Df,Dp,A,B,C}(cell_trian,face_trian,glue)
  end
end

get_background_model(t::HybridTriangulation) = get_background_model(t.cell_trian)
get_grid(t::HybridTriangulation) = get_grid(t.face_trian)
get_glue(t::HybridTriangulation,::Val{d}) where d = get_glue(t.face_trian,Val(d))
get_facet_normal(t::HybridTriangulation) = get_facet_normal(t.face_trian)

function HybridTriangulation(model::DiscreteModel{Dc},cell_to_bgcell) where Dc
  Df = Dc-1 # TODO: Expand to arbitrary face dimension
  glue = CellToFaceGlue(model,cell_to_bgcell,Df)

  face_to_bgface = glue.face_glue.face_to_bgface
  cell_trian = Triangulation(model,cell_to_bgcell)
  face_trian = view(Triangulation(ReferenceFE{Df},model),face_to_bgface)
  
  return HybridTriangulation(cell_trian,face_trian,glue)
end
