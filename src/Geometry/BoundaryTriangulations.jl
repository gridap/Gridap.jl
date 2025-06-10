
struct FaceToCellGlue{A,B,C,D} <: GridapType
  face_to_bgface::A
  bgface_to_lcell::B
  face_to_cell::Vector{Int32}
  face_to_lface::Vector{Int8}
  face_to_lcell::Vector{Int8}
  face_to_ftype::C
  cell_to_ctype::D
  cell_to_lface_to_pindex::Table{Int8,Vector{Int8},Vector{Int32}}
  ctype_to_lface_to_ftype::Vector{Vector{Int8}}
end

function FaceToCellGlue(
  topo::GridTopology,
  cell_grid::Grid,
  face_grid::Grid,
  face_to_bgface::AbstractVector,
  bgface_to_lcell::AbstractVector)

  D = num_cell_dims(cell_grid)
  cD = num_cell_dims(face_grid)
  bgface_to_cells = get_faces(topo,cD,D)
  cell_to_bgfaces = get_faces(topo,D,cD)
  cell_to_lface_to_pindex = Table(get_cell_permutations(topo,cD))

  bgface_to_cell = lazy_map(getindex,bgface_to_cells, bgface_to_lcell)
  bgface_to_lface = find_local_index(bgface_to_cell, cell_to_bgfaces)

  face_to_cell = collect(Int32,lazy_map(Reindex(bgface_to_cell), face_to_bgface))
  face_to_lface = collect(Int8,lazy_map(Reindex(bgface_to_lface), face_to_bgface))
  face_to_lcell = collect(Int8,lazy_map(Reindex(bgface_to_lcell), face_to_bgface))

  f = (p)->fill(Int8(UNSET),num_faces(p,cD))
  ctype_to_lface_to_ftype = map( f, get_reffes(cell_grid) )
  face_to_ftype = get_cell_type(face_grid)
  cell_to_ctype = get_cell_type(cell_grid)

  _fill_ctype_to_lface_to_ftype!(
    ctype_to_lface_to_ftype,
    face_to_cell,
    face_to_lface,
    face_to_ftype,
    cell_to_ctype)

  FaceToCellGlue(
    face_to_bgface,
    bgface_to_lcell,
    face_to_cell,
    face_to_lface,
    face_to_lcell,
    face_to_ftype,
    cell_to_ctype,
    cell_to_lface_to_pindex,
    ctype_to_lface_to_ftype)
end

function _fill_ctype_to_lface_to_ftype!(
  ctype_to_lface_to_ftype,
  face_to_cell,
  face_to_lface,
  face_to_ftype,
  cell_to_ctype)
  for (face, cell) in enumerate(face_to_cell)
    ctype = cell_to_ctype[cell]
    lface = face_to_lface[face]
    ftype = face_to_ftype[face]
    ctype_to_lface_to_ftype[ctype][lface] = ftype
  end
end

function get_children(n::TreeNode, a::FaceToCellGlue)
  (
    similar_tree_node(n,a.face_to_cell),
    similar_tree_node(n,a.face_to_lface),
    similar_tree_node(n,a.face_to_lcell),
    similar_tree_node(n,a.face_to_ftype),
    similar_tree_node(n,a.cell_to_ctype),
    similar_tree_node(n,a.cell_to_lface_to_pindex),
    similar_tree_node(n,a.ctype_to_lface_to_ftype)
  )
end

"""
"""
struct BoundaryTriangulation{Dc,Dp,A,B} <: Triangulation{Dc,Dp}
  trian::A
  glue::B
  function BoundaryTriangulation(
    trian::BodyFittedTriangulation,
    glue)

    Dc = num_cell_dims(trian)
    Dp = num_point_dims(trian)
    A = typeof(trian)
    B = typeof(glue)
    new{Dc,Dp,A,B}(trian,glue)
  end
end

"""
    Boundary(args...; kwargs...)

Alias for [`BoundaryTriangulation`](@ref)(args...; kwargs...).
"""
function Boundary(args...;kwargs...)
  BoundaryTriangulation(args...;kwargs...)
end

# Constructors

"""
    BoundaryTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
    BoundaryTriangulation(model::DiscreteModel)
"""
function BoundaryTriangulation(
  model::DiscreteModel,
  face_to_bgface::AbstractVector{<:Integer},
  bgface_to_lcell::AbstractVector{<:Integer})

  D = num_cell_dims(model)
  topo = get_grid_topology(model)
  bgface_grid = Grid(ReferenceFE{D-1},model)

  face_grid = view(bgface_grid,face_to_bgface)
  cell_grid = get_grid(model)
  glue = FaceToCellGlue(topo,cell_grid,face_grid,face_to_bgface,bgface_to_lcell)
  trian = BodyFittedTriangulation(model,face_grid,face_to_bgface)

  BoundaryTriangulation(trian,glue)
end

function BoundaryTriangulation(
  model::DiscreteModel,
  face_to_bgface::AbstractVector{<:Integer})
  BoundaryTriangulation(model,face_to_bgface,Fill(1,num_facets(model)))
end

function BoundaryTriangulation(
  model::DiscreteModel,
  bgface_to_mask::AbstractVector{Bool},
  bgface_to_lcell::AbstractVector{<:Integer})

  face_to_bgface = findall(bgface_to_mask)
  BoundaryTriangulation(model,face_to_bgface,bgface_to_lcell)
end

function BoundaryTriangulation(
  model::DiscreteModel, bgface_to_mask::AbstractVector{Bool}, lcell::Integer=1)
  BoundaryTriangulation(model,bgface_to_mask, Fill(lcell,num_facets(model)) )
end

"""
    BoundaryTriangulation(model::DiscreteModel,labeling::FaceLabeling;tags::Vector{Int})
    BoundaryTriangulation(model::DiscreteModel,labeling::FaceLabeling;tags::Vector{String})
    BoundaryTriangulation(model::DiscreteModel,labeling::FaceLabeling;tag::Int)
    BoundaryTriangulation(model::DiscreteModel,labeling::FaceLabeling;tag::String)
"""
function BoundaryTriangulation(model::DiscreteModel,labeling::FaceLabeling;tags=nothing)
  D = num_cell_dims(model)
  if tags == nothing
    topo = get_grid_topology(model)
    face_to_mask = get_isboundary_face(topo,D-1)
  else
    face_to_mask = get_face_mask(labeling,tags,D-1)
  end
  BoundaryTriangulation(model,face_to_mask)
end

"""
    BoundaryTriangulation(model::DiscreteModel,tags::Vector{Int})
    BoundaryTriangulation(model::DiscreteModel,tags::Vector{String})
    BoundaryTriangulation(model::DiscreteModel,tag::Int)
    BoundaryTriangulation(model::DiscreteModel,tag::String)
"""
function BoundaryTriangulation(model::DiscreteModel;tags=nothing)
  labeling = get_face_labeling(model)
  BoundaryTriangulation(model,labeling,tags=tags)
end

function BoundaryTriangulation(rtrian::Triangulation,args...;kwargs...)
  rmodel = get_active_model(rtrian)
  dtrian = BoundaryTriangulation(rmodel,args...;kwargs...)
  CompositeTriangulation(rtrian,dtrian)
end

# API

get_background_model(t::BoundaryTriangulation) = get_background_model(t.trian)
get_grid(t::BoundaryTriangulation) = get_grid(t.trian)
get_glue(t::BoundaryTriangulation{D},::Val{D}) where D = get_glue(t.trian,Val(D))

function get_glue(trian::BoundaryTriangulation,::Val{Dp}) where Dp
  model = get_background_model(trian)
  Dm = num_cell_dims(model)
  get_glue(trian,Val(Dp),Val(Dm))
end

function get_glue(trian::BoundaryTriangulation,::Val{Dp},::Val{Dm}) where {Dp,Dm}
  nothing
end

function get_glue(trian::BoundaryTriangulation,::Val{D},::Val{D}) where D
  tface_to_mface = trian.glue.face_to_cell
  face_to_q_vertex_coords = _compute_face_to_q_vertex_coords(trian)
  f(p) = get_shapefuns(LagrangianRefFE(Float64,get_polytope(p),1))
  ftype_to_shapefuns = map( f, get_reffes(trian) )
  face_to_shapefuns = expand_cell_data(ftype_to_shapefuns,trian.glue.face_to_ftype)
  face_s_q = lazy_map(linear_combination,face_to_q_vertex_coords,face_to_shapefuns)
  tface_to_mface_map = face_s_q
  mface_to_tface = nothing
  FaceToFaceGlue(tface_to_mface,tface_to_mface_map,mface_to_tface)
end

function get_facet_normal(trian::BoundaryTriangulation, boundary_trian_glue::FaceToCellGlue)
  cell_grid = get_grid(get_background_model(trian.trian))

  ## Reference normal
  function f(r)
    p = get_polytope(r)
    lface_to_n = get_facet_normal(p)
    lface_to_pindex_to_perm = get_face_vertex_permutations(p,num_cell_dims(p)-1)
    nlfaces = length(lface_to_n)
    lface_pindex_to_n = [ fill(lface_to_n[lface],length(lface_to_pindex_to_perm[lface])) for lface in 1:nlfaces ]
    lface_pindex_to_n
  end
  ctype_lface_pindex_to_nref = map(f, get_reffes(cell_grid))
  face_to_nref = FaceCompressedVector(ctype_lface_pindex_to_nref,boundary_trian_glue)
  face_s_nref = lazy_map(constant_field,face_to_nref)

  # Inverse of the Jacobian transpose
  cell_q_x = get_cell_map(cell_grid)
  cell_q_Jt = lazy_map(∇,cell_q_x)
  cell_q_invJt = lazy_map(Operation(pinvJt),cell_q_Jt)
  face_q_invJt = lazy_map(Reindex(cell_q_invJt),boundary_trian_glue.face_to_cell)

  # Change of domain
  D = num_cell_dims(cell_grid)
  boundary_trian_glue = get_glue(trian,Val(D))
  face_s_q = boundary_trian_glue.tface_to_mface_map
  face_s_invJt = lazy_map(∘,face_q_invJt,face_s_q)
  face_s_n = lazy_map(Broadcasting(Operation(push_normal)),face_s_invJt,face_s_nref)
  Fields.MemoArray(face_s_n)
end

function get_edge_tangent(trian::BoundaryTriangulation, boundary_trian_glue::FaceToCellGlue)

  face_s_n = get_facet_normal(trian, boundary_trian_glue)

  rotate_normal_2d(n) = VectorValue(n[2], -n[1])

  face_s_t = lazy_map(Operation(rotate_normal_2d), face_s_n)

  return Fields.MemoArray(face_s_t)
end

function get_edge_tangent(trian::BoundaryTriangulation{Dc,Dp,A,<:FaceToCellGlue}) where {Dc,Dp,A}
  Dp != 2 && error("get_edge_tangent only implemented for 2D")
  glue = trian.glue
  get_edge_tangent(trian, glue)
end

function get_facet_normal(trian::BoundaryTriangulation{Dc,Dp,A,<:FaceToCellGlue}) where {Dc,Dp,A}
  glue = trian.glue
  get_facet_normal(trian, glue)
end

function push_normal(invJt,n)
  v = invJt⋅n
  m = sqrt(inner(v,v))
  if m < eps()
    return zero(v)
  else
    return v/m
  end
end

function _compute_face_to_q_vertex_coords(trian::BoundaryTriangulation,glue)
  d = num_cell_dims(trian)
  cell_grid = get_grid(get_background_model(trian.trian))
  polytopes = map(get_polytope, get_reffes(cell_grid))
  cell_to_ctype = glue.cell_to_ctype
  ctype_to_lvertex_to_qcoords = map(get_vertex_coordinates, polytopes)
  ctype_to_lface_to_lvertices = map((p)->get_faces(p,d,0), polytopes)
  ctype_to_lface_to_pindex_to_perm = map( (p)->get_face_vertex_permutations(p,d), polytopes)

  P = eltype(eltype(ctype_to_lvertex_to_qcoords))
  D = num_components(P)
  T = eltype(P)
  ctype_to_lface_to_pindex_to_qcoords = Vector{Vector{Vector{Point{D,T}}}}[]

  for (ctype, lface_to_pindex_to_perm) in enumerate(ctype_to_lface_to_pindex_to_perm)
    lvertex_to_qcoods = ctype_to_lvertex_to_qcoords[ctype]
    lface_to_pindex_to_qcoords = Vector{Vector{Point{D,T}}}[]
    for (lface, pindex_to_perm) in enumerate(lface_to_pindex_to_perm)
      cfvertex_to_lvertex = ctype_to_lface_to_lvertices[ctype][lface]
      nfvertices = length(cfvertex_to_lvertex)
      pindex_to_qcoords = Vector{Vector{Point{D,T}}}(undef,length(pindex_to_perm))
      for (pindex, cfvertex_to_ffvertex) in enumerate(pindex_to_perm)
        ffvertex_to_qcoords = zeros(Point{D,T},nfvertices)
        for (cfvertex, ffvertex) in enumerate(cfvertex_to_ffvertex)
          lvertex = cfvertex_to_lvertex[cfvertex]
          qcoords = lvertex_to_qcoods[lvertex]
          ffvertex_to_qcoords[ffvertex] = qcoords
        end
        pindex_to_qcoords[pindex] = ffvertex_to_qcoords
      end
      push!(lface_to_pindex_to_qcoords,pindex_to_qcoords)
    end
    push!(ctype_to_lface_to_pindex_to_qcoords,lface_to_pindex_to_qcoords)
  end

  FaceCompressedVector(ctype_to_lface_to_pindex_to_qcoords,glue)
end

function _compute_face_to_q_vertex_coords(trian::BoundaryTriangulation{Dc,Dp,A,<:FaceToCellGlue}) where {Dc,Dp,A}
  _compute_face_to_q_vertex_coords(trian,trian.glue)
end

struct FaceCompressedVector{T,G<:FaceToCellGlue} <: AbstractVector{T}
  ctype_lface_pindex_to_value::Vector{Vector{Vector{T}}}
  glue::G
end

Base.size(a::FaceCompressedVector) = (length(a.glue.face_to_cell),)

Base.IndexStyle(::Type{<:FaceCompressedVector}) = IndexLinear()

function Base.getindex(a::FaceCompressedVector,face::Integer)
  cell = a.glue.face_to_cell[face]
  ctype = a.glue.cell_to_ctype[cell]
  lface = a.glue.face_to_lface[face]
  p = a.glue.cell_to_lface_to_pindex.ptrs[cell]-1
  pindex = a.glue.cell_to_lface_to_pindex.data[p+lface]
  value = a.ctype_lface_pindex_to_value[ctype][lface][pindex]
  value
end

function get_children(n::TreeNode, a::FaceCompressedVector)
  (similar_tree_node(n,a.ctype_lface_pindex_to_value),similar_tree_node(n,a.glue))
end

function lazy_map(k::Fields.LinearCombinationMap,::Type{T},b::FaceCompressedVector,c::Fill) where T
  d = CompressedArray([c.value,],Fill(1,length(c)))
  lazy_map(k,T,a,b,d)
end

function lazy_map(k::Fields.LinearCombinationMap,::Type{T},b::FaceCompressedVector,c::CompressedArray) where T
  if c.ptrs === b.glue.face_to_ftype || c.ptrs == b.glue.face_to_ftype

    ctype_lface_pindex_to_r = Vector{Vector{Vector{T}}}(undef,length(b.ctype_lface_pindex_to_value))
    for (ctype, lface_pindex_to_value) in enumerate(b.ctype_lface_pindex_to_value)
      lface_pindex_to_r = Vector{Vector{T}}(undef,length(lface_pindex_to_value))
      for (lface, pindex_to_value) in enumerate(lface_pindex_to_value)
        pindex_to_r = Vector{T}(undef,length(pindex_to_value))
        ftype = b.glue.ctype_to_lface_to_ftype[ctype][lface]
        if ftype != UNSET
          for (pindex, value) in enumerate(pindex_to_value)
            pindex_to_r[pindex] = k(value,c.values[ftype])
          end
        end
        lface_pindex_to_r[lface] = pindex_to_r
      end
      ctype_lface_pindex_to_r[ctype] = lface_pindex_to_r
    end
    FaceCompressedVector(ctype_lface_pindex_to_r,b.glue)

  else
    @notimplemented
  end
end

function lazy_map(k::typeof(evaluate),::Type{T},a::Fill,b::FaceCompressedVector) where T
  @check length(a) == length(b)
  ctype_lface_pindex_to_r = Vector{Vector{Vector{T}}}(undef,length(b.ctype_lface_pindex_to_value))
  for (ctype, lface_pindex_to_value) in enumerate(b.ctype_lface_pindex_to_value)
    lface_pindex_to_r = Vector{Vector{T}}(undef,length(lface_pindex_to_value))
    for (lface, pindex_to_value) in enumerate(lface_pindex_to_value)
      pindex_to_r = Vector{T}(undef,length(pindex_to_value))
      ftype = b.glue.ctype_to_lface_to_ftype[ctype][lface]
      if ftype != UNSET
        for (pindex, value) in enumerate(pindex_to_value)
          pindex_to_r[pindex] = evaluate(a.value,value)
        end
      end
      lface_pindex_to_r[lface] = pindex_to_r
    end
    ctype_lface_pindex_to_r[ctype] = lface_pindex_to_r
  end
  FaceCompressedVector(ctype_lface_pindex_to_r,b.glue)
end

function lazy_map(
  k::typeof(evaluate),
  ::Type{T},
  a::CompressedArray,
  b::FaceCompressedVector) where T

  @check length(a) == length(b)
  ctype_lface_pindex_to_r = Vector{Vector{Vector{T}}}(undef,length(b.ctype_lface_pindex_to_value))
  for (ctype, lface_pindex_to_value) in enumerate(b.ctype_lface_pindex_to_value)
    lface_pindex_to_r = Vector{Vector{T}}(undef,length(lface_pindex_to_value))
    for (lface, pindex_to_value) in enumerate(lface_pindex_to_value)
      pindex_to_r = Vector{T}(undef,length(pindex_to_value))
      ftype = b.glue.ctype_to_lface_to_ftype[ctype][lface]
      if ftype != UNSET
        for (pindex, value) in enumerate(pindex_to_value)
          pindex_to_r[pindex] = evaluate(a.values[ctype],value)
        end
      end
      lface_pindex_to_r[lface] = pindex_to_r
    end
    ctype_lface_pindex_to_r[ctype] = lface_pindex_to_r
  end
  FaceCompressedVector(ctype_lface_pindex_to_r,b.glue)
end

function lazy_map(k::typeof(evaluate),::Type{T},b::FaceCompressedVector,a::Fill) where T
  @check length(a) == length(b)
  ctype_lface_pindex_to_r = Vector{Vector{Vector{T}}}(undef,length(b.ctype_lface_pindex_to_value))
  for (ctype, lface_pindex_to_value) in enumerate(b.ctype_lface_pindex_to_value)
    lface_pindex_to_r = Vector{Vector{T}}(undef,length(lface_pindex_to_value))
    for (lface, pindex_to_value) in enumerate(lface_pindex_to_value)
      pindex_to_r = Vector{T}(undef,length(pindex_to_value))
      ftype = b.glue.ctype_to_lface_to_ftype[ctype][lface]
      if ftype != UNSET
        for (pindex, value) in enumerate(pindex_to_value)
          pindex_to_r[pindex] = evaluate(value,a.value)
        end
      end
      lface_pindex_to_r[lface] = pindex_to_r
    end
    ctype_lface_pindex_to_r[ctype] = lface_pindex_to_r
  end
  FaceCompressedVector(ctype_lface_pindex_to_r,b.glue)
end

function lazy_map(k::typeof(evaluate),::Type{T},a::Fill,b::FaceCompressedVector,c::FaceCompressedVector) where T
  if b.glue !== c.glue
    return LazyArray(T,a,b,c)
  end
  @check length(a) == length(b)
  ctype_lface_pindex_to_r = Vector{Vector{Vector{T}}}(undef,length(b.ctype_lface_pindex_to_value))
  for (ctype, lface_pindex_to_value) in enumerate(b.ctype_lface_pindex_to_value)
    lface_pindex_to_r = Vector{Vector{T}}(undef,length(lface_pindex_to_value))
    for (lface, pindex_to_value) in enumerate(lface_pindex_to_value)
      pindex_to_r = Vector{T}(undef,length(pindex_to_value))
      ftype = b.glue.ctype_to_lface_to_ftype[ctype][lface]
      if ftype != UNSET
        for (pindex, bvalue) in enumerate(pindex_to_value)
          cvalue = c.ctype_lface_pindex_to_value[ctype][lface][pindex]
          pindex_to_r[pindex] = evaluate(a.value,bvalue,cvalue)
        end
      end
      lface_pindex_to_r[lface] = pindex_to_r
    end
    ctype_lface_pindex_to_r[ctype] = lface_pindex_to_r
  end
  FaceCompressedVector(ctype_lface_pindex_to_r,b.glue)
end

