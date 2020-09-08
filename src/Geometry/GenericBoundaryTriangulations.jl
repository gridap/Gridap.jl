
struct FaceToCellGlue{O} <: GridapType
  orientation::Val{O}
  face_to_cell::Vector{Int}
  face_to_lface::Vector{Int8}
  cell_to_ctype::Vector{Int8}
  cell_to_lface_to_pindex::Table{Int8,Vector{Int8},Vector{Int32}}
  ctype_to_lface_to_ftype::Vector{Vector{Int8}}
end

is_oriented(::FaceToCellGlue{O}) where O = O

"""
    struct GenericBoundaryTriangulation{Dc,Dp,Gf,Gc,O} <: BoundaryTriangulation{Dc,Dp}
      face_trian::Gf
      cell_trian::Gc
      # + private fields
    end
"""
struct GenericBoundaryTriangulation{Dc,Dp,Gf,Gc,O,A} <: BoundaryTriangulation{Dc,Dp}
  face_trian::Gf
  cell_trian::Gc
  glue::FaceToCellGlue{O}
  cell_around::A
  memo::Dict

  function GenericBoundaryTriangulation(
    face_trian::Triangulation,
    cell_trian::Triangulation,
    glue::FaceToCellGlue{O},
    cell_around=1) where O

    Dc = num_cell_dims(face_trian)
    Dp = num_point_dims(face_trian)
    Gf = typeof(face_trian)
    Gc = typeof(cell_trian)
    A = typeof(cell_around)
    new{Dc,Dp,Gf,Gc,O,A}(face_trian, cell_trian, glue, cell_around, Dict())
  end

end

get_memo(a::GenericBoundaryTriangulation) = a.memo

function GenericBoundaryTriangulation(
  face_trian::Triangulation,
  cell_trian::Triangulation,
  topo::GridTopology,
  face_to_oldface::Vector{Int},
  icell_around=1)


  D = num_cell_dims(cell_trian)
  oldface_to_cells = get_faces(topo,D-1,D)
  cell_to_oldfaces = get_faces(topo,D,D-1)
  cell_to_lface_to_pindex = Table(get_cell_permutations(topo,D-1))
  orientation = OrientationStyle(topo)

  oldface_to_cell = get_local_item(oldface_to_cells, icell_around)
  oldface_to_lface = find_local_index(oldface_to_cell, cell_to_oldfaces)
  face_to_cell = collect(Int,reindex(oldface_to_cell, face_to_oldface))
  face_to_lface = collect(Int8,reindex(oldface_to_lface, face_to_oldface))

  f = (p)->fill(Int8(UNSET),num_faces(p,D-1))
  ctype_to_lface_to_ftype = map( f, get_reffes(cell_trian) )
  face_to_ftype = get_cell_type(face_trian)
  cell_to_ctype = collect(Int8,get_cell_type(cell_trian))
  _fill_ctype_to_lface_to_ftype!(
    ctype_to_lface_to_ftype,
    face_to_cell,
    face_to_lface,
    face_to_ftype,
    cell_to_ctype)

  glue = FaceToCellGlue(
    orientation,
    face_to_cell,
    face_to_lface,
    cell_to_ctype,
    cell_to_lface_to_pindex,
    ctype_to_lface_to_ftype)

  GenericBoundaryTriangulation(face_trian, cell_trian, glue, icell_around)

end

"""
    GenericBoundaryTriangulation(model::DiscreteModel,face_to_mask::Vector{Bool})
"""
function GenericBoundaryTriangulation(
  model::DiscreteModel,face_to_mask::Vector{Bool},icell_around=1)
  D = num_cell_dims(model)
  topo = get_grid_topology(model)
  oldface_grid = Grid(ReferenceFE{D-1},model)
  cell_grid = Grid(ReferenceFE{D},model)
  face_to_oldface = findall(face_to_mask)
  face_trian = TriangulationPortion(oldface_grid,face_to_oldface)
  GenericBoundaryTriangulation(face_trian,cell_grid,topo,face_to_oldface,icell_around)
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

# Triangulation Interface

function get_cell_coordinates(trian::GenericBoundaryTriangulation)
  get_cell_coordinates(trian.face_trian)
end

function get_reffes(trian::GenericBoundaryTriangulation)
    get_reffes(trian.face_trian)
end

function get_cell_type(trian::GenericBoundaryTriangulation)
  get_cell_type(trian.face_trian)
end

function get_face_to_cell(trian::GenericBoundaryTriangulation)
  trian.glue.face_to_cell
end

function get_face_to_lface(trian::GenericBoundaryTriangulation)
  trian.glue.face_to_lface
end

function get_volume_triangulation(trian::GenericBoundaryTriangulation)
  trian.cell_trian
end

function get_face_to_face(trian::GenericBoundaryTriangulation)
  trian.face_trian.cell_to_oldcell
end

function get_cell_around(trian::GenericBoundaryTriangulation)
  trian.cell_around
end

# Cell to face map

function compute_face_to_cell_map(trian::GenericBoundaryTriangulation)
  face_to_fvertex_to_qcoors = FaceCellCoordinates(trian)
  f = (p)-> get_shapefuns(LagrangianRefFE(Float64,get_polytope(p),1))
  ftype_to_shapefuns = map( f, get_reffes(trian.face_trian) )
  face_to_ftype = collect(Int8, get_cell_type(trian.face_trian))
  face_to_shapefuns = CompressedArray(ftype_to_shapefuns, face_to_ftype)
  array = lincomb(face_to_shapefuns,face_to_fvertex_to_qcoors)
  GenericCellField(array)
end

struct FaceCellCoordinates{D,T,O} <: AbstractVector{Vector{Point{D,T}}}
  glue::FaceToCellGlue{O}
  cell_to_ctype::Vector{Int8}
  ctype_to_lvertex_to_qcoords::Vector{Vector{Point{D,T}}}
  ctype_to_lface_to_lvertices::Vector{Vector{Vector{Int}}}
  ctype_to_lface_to_pindex_to_perm::Vector{Vector{Vector{Vector{Int}}}}
  ctype_to_lface_to_pindex_to_qcoords::Vector{Vector{Vector{Vector{Point{D,T}}}}}
  function FaceCellCoordinates(trian::GenericBoundaryTriangulation)
  
    d = num_cell_dims(trian)
    polytopes = map(get_polytope, get_reffes(trian.cell_trian))
    cell_to_ctype = trian.glue.cell_to_ctype
    ctype_to_lvertex_to_qcoords = map(get_vertex_coordinates, polytopes)
    ctype_to_lface_to_lvertices = map((p)->get_faces(p,d,0), polytopes)
    ctype_to_lface_to_pindex_to_perm = map( (p)->get_face_vertex_permutations(p,d), polytopes)
  
    P = eltype(eltype(ctype_to_lvertex_to_qcoords))
    D = length(P)
    T = eltype(P)
    O = is_oriented(trian.glue)
  
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

    new{D,T,O}(
      trian.glue,
      cell_to_ctype,
      ctype_to_lvertex_to_qcoords,
      ctype_to_lface_to_lvertices,
      ctype_to_lface_to_pindex_to_perm,
      ctype_to_lface_to_pindex_to_qcoords)
  
  end
end

Base.IndexStyle(::Type{FaceCellCoordinates{D,T,O}}) where {D,T,O} = IndexLinear()

Base.size(a::FaceCellCoordinates) = (length(a.glue.face_to_cell),)

function Base.getindex(a::FaceCellCoordinates,face::Integer)
  cell = a.glue.face_to_cell[face]
  lface = a.glue.face_to_lface[face]
  ctype = a.cell_to_ctype[cell]
  p = a.glue.cell_to_lface_to_pindex.ptrs[cell]-1
  pindex = a.glue.cell_to_lface_to_pindex.data[p+lface]
  a.ctype_to_lface_to_pindex_to_qcoords[ctype][lface][pindex]
end

# Value of face to cell map

function apply_lincomb(ax::CompressedArray,b::FaceCellCoordinates)
  FaceToCellMapValue(ax,b)
end

struct FaceToCellMapValue{D,T,O} <: AbstractVector{Vector{Point{D,T}}}
  b::FaceCellCoordinates{D,T,O}
  ctype_to_lface_to_pindex_to_qpoints::Vector{Vector{Vector{Vector{Point{D,T}}}}}

  function FaceToCellMapValue(ax::CompressedArray,b::FaceCellCoordinates{D,T,O}) where {D,T,O}

    # Precompute the integration point coordinates for all possible cases

    ctype_to_lface_to_ftype = b.glue.ctype_to_lface_to_ftype
    ctype_to_lface_to_lvertices = b.ctype_to_lface_to_lvertices
    ctype_to_lvertex_to_qcoords = b.ctype_to_lvertex_to_qcoords
    ctype_to_lface_to_pindex_to_perm = b.ctype_to_lface_to_pindex_to_perm
    ftype_to_shapefunsvals = ax.values

    # Allocate
    ctype_to_lface_to_pindex_to_qpoints = Vector{Vector{Vector{Point{D,T}}}}[]
    for (ctype, lface_to_pindex_to_perm) in enumerate(ctype_to_lface_to_pindex_to_perm)
      lface_to_pindex_to_qpoints = Vector{Vector{Point{D,T}}}[]
      for (lface, pindex_to_perm) in enumerate(lface_to_pindex_to_perm)
        pindex_to_qpoints = Vector{Vector{Point{D,T}}}(undef,length(pindex_to_perm))
        push!(lface_to_pindex_to_qpoints,pindex_to_qpoints)
      end
      push!(ctype_to_lface_to_pindex_to_qpoints,lface_to_pindex_to_qpoints)
    end

    # Fill
    for (ctype, lface_to_ftype) in enumerate(ctype_to_lface_to_ftype)
      for (lface, ftype) in enumerate(lface_to_ftype)
        if ftype == UNSET
          continue
        end

        cfvertex_to_lvertex = ctype_to_lface_to_lvertices[ctype][lface]
        pindex_to_perm = ctype_to_lface_to_pindex_to_perm[ctype][lface]
        lvertex_to_qcoods = ctype_to_lvertex_to_qcoords[ctype]
        nfvertices = length(cfvertex_to_lvertex)
        ffvertex_to_qcoords = zeros(Point{D,T},nfvertices)
        shapefunsvals = ftype_to_shapefunsvals[ftype]

        for (pindex, cfvertex_to_ffvertex) in enumerate(pindex_to_perm)
          for (cfvertex, ffvertex) in enumerate(cfvertex_to_ffvertex)
            lvertex = cfvertex_to_lvertex[cfvertex]
            qcoords = lvertex_to_qcoods[lvertex]
            ffvertex_to_qcoords[ffvertex] = qcoords
          end
          qpoints = shapefunsvals*ffvertex_to_qcoords
          ctype_to_lface_to_pindex_to_qpoints[ctype][lface][pindex] = qpoints
        end
      end
    end

    new{D,T,O}(b,ctype_to_lface_to_pindex_to_qpoints)

  end
end

Base.size(v::FaceToCellMapValue) = (length(v.b.glue.face_to_cell),)

Base.IndexStyle(::Type{FaceToCellMapValue{C,T,O}}) where {C,T,O} = IndexLinear()

function Base.getindex(a::FaceToCellMapValue,face::Integer)
  cell = a.b.glue.face_to_cell[face]
  lface = a.b.glue.face_to_lface[face]
  ctype = a.b.cell_to_ctype[cell]
  p = a.b.glue.cell_to_lface_to_pindex.ptrs[cell]-1
  pindex = a.b.glue.cell_to_lface_to_pindex.data[p+lface]
  a.ctype_to_lface_to_pindex_to_qpoints[ctype][lface][pindex]
end

function evaluate_field_array(g::CompressedArray{<:Field},fx::FaceToCellMapValue)
  CellShapeFunsAtFaces(g,fx)
end

struct CellShapeFunsAtFaces{D,T,O,V} <: AbstractVector{V}
  b::FaceCellCoordinates{D,T,O}
  ctype_to_lface_to_pindex_to_shfnsvals::Vector{Vector{Vector{V}}}
  function CellShapeFunsAtFaces(g::CompressedArray{<:Field},fx::FaceToCellMapValue{D,T,O}) where {D,T,O}

    # Precompute shapefuns for all cases

    ctype_to_shapefuns = g.values
    face_to_ctype = g.ptrs
    ctype_to_lface_to_pindex_to_perm = fx.b.ctype_to_lface_to_pindex_to_perm
    ctype_to_lface_to_pindex_to_qpoints = fx.ctype_to_lface_to_pindex_to_qpoints
    ctype_to_lface_to_ftype = fx.b.glue.ctype_to_lface_to_ftype

    # Allocate
    sfuns = first(ctype_to_shapefuns)
    V = field_return_type(sfuns,zeros(Point{D,T},1))
    ctype_to_lface_to_pindex_to_shfnsvals = Vector{Vector{V}}[]
    for (ctype, lface_to_pindex_to_perm) in enumerate(ctype_to_lface_to_pindex_to_perm)
      lface_to_pindex_to_shfnsvals = Vector{V}[]
      for (lface, pindex_to_perm) in enumerate(lface_to_pindex_to_perm)
        pindex_to_shfnsvals = Vector{V}(undef,length(pindex_to_perm))
        push!(lface_to_pindex_to_shfnsvals,pindex_to_shfnsvals)
      end
      push!(ctype_to_lface_to_pindex_to_shfnsvals,lface_to_pindex_to_shfnsvals)
    end

    # Fill
    for (ctype, lface_to_ftype) in enumerate(ctype_to_lface_to_ftype)
      for (lface, ftype) in enumerate(lface_to_ftype)
        if ftype == UNSET
          continue
        end
        pindex_to_perm = ctype_to_lface_to_pindex_to_perm[ctype][lface]
        sfuns = ctype_to_shapefuns[ctype]
        for (pindex, cfvertex_to_ffvertex) in enumerate(pindex_to_perm)
          qpoints = ctype_to_lface_to_pindex_to_qpoints[ctype][lface][pindex]
          shfnsvals = evaluate(sfuns,qpoints)
          ctype_to_lface_to_pindex_to_shfnsvals[ctype][lface][pindex] = shfnsvals
        end
      end
    end

    new{D,T,O,V}(fx.b,ctype_to_lface_to_pindex_to_shfnsvals)

  end
end

Base.size(v::CellShapeFunsAtFaces) = (length(v.b.glue.face_to_cell),)

Base.IndexStyle(::Type{CellShapeFunsAtFaces{C,T,O,V}}) where {C,T,O,V} = IndexLinear()

function Base.getindex(a::CellShapeFunsAtFaces,face::Integer)
  cell = a.b.glue.face_to_cell[face]
  lface = a.b.glue.face_to_lface[face]
  ctype = a.b.cell_to_ctype[cell]
  p = a.b.glue.cell_to_lface_to_pindex.ptrs[cell]-1
  pindex = a.b.glue.cell_to_lface_to_pindex.data[p+lface]
  a.ctype_to_lface_to_pindex_to_shfnsvals[ctype][lface][pindex]
end

# Normal vector

function get_normal_vector(trian::GenericBoundaryTriangulation)
  k = NormalVectorValued()
  refn = ReferenceNormal(trian)
  ϕv = get_cell_map(trian.cell_trian)
  ϕ = get_cell_map(trian)
  J = (∇(ϕv)∘inverse_map(ϕv))∘ϕ
  a = apply(k,get_array(J),refn)
  ϕ = get_cell_map(trian)
  GenericCellField(a)∘inverse_map(ϕ)
end

struct ReferenceNormal{D,T} <: AbstractVector{Point{D,T}}
  face_to_cell::Vector{Int}
  face_to_lface::Vector{Int8}
  cell_to_ctype::Vector{Int8}
  ctype_to_lface_to_qnormal::Vector{Vector{Point{D,T}}}
end

function ReferenceNormal(trian::GenericBoundaryTriangulation)
  face_to_cell = trian.glue.face_to_cell
  face_to_lface = trian.glue.face_to_lface
  cell_to_ctype = trian.glue.cell_to_ctype
  f = (r) -> get_facet_normals(get_polytope(r))
  ctype_to_lface_to_qnormal = map(f, get_reffes(trian.cell_trian))
  ReferenceNormal(
    face_to_cell,
    face_to_lface,
    cell_to_ctype,
    ctype_to_lface_to_qnormal)
end

Base.size(a::ReferenceNormal) = (length(a.face_to_cell),)

Base.IndexStyle(::Type{ReferenceNormal{D,T}}) where {D,T} = IndexLinear()

function Base.getindex(a::ReferenceNormal,face::Integer)
  cell = a.face_to_cell[face]
  lface = a.face_to_lface[face]
  ctype = a.cell_to_ctype[cell]
  a.ctype_to_lface_to_qnormal[ctype][lface]
end

struct NormalField <: Field end

struct NormalVectorValued <: Kernel end

function apply_kernel!(cache,k::NormalVectorValued,J,nref)
  NormalField()
end

function kernel_evaluate(k::NormalVectorValued,x,J,refn)
  Jx = evaluate_field_array(J,x)
  k = bcast(_map_normal)
  apply(k,Jx,refn)
end

function _map_normal(J::TensorValue{D,D,T},n::VectorValue{D,T}) where {D,T}
  v = inv(J)⋅n
  m = sqrt(inner(v,v))
  if m < eps()
    return zero(n)
  else
    return v/m
  end
end

