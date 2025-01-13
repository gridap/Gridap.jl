
"""
    struct PatchTopology{Dc,Dp} <: GridapType
      topo :: GridTopology{Dc,Dp}
      d_to_patch_to_dpfaces :: Vector{Table{Int32,Vector{Int32},Vector{Int32}}}
      d_to_dpface_to_dface :: Vector{Vector{Int32}}
    end

# Fields: 

- `topo`: Underlying grid topology
- `d_to_patch_to_dpfaces`: For each dimension `d`, a table mapping each patch to it's d-faces.
"""
struct PatchTopology{Dc,Dp} <: GridapType
  topo :: GridTopology{Dc,Dp}
  d_to_patch_to_dfaces :: Vector{Table{Int32,Vector{Int32},Vector{Int32}}}
end

function PatchTopology(topo::GridTopology{Dc},patch_cells::Table) where Dc
  d_to_patch_to_dpfaces = Vector{Table{Int32,Vector{Int32},Vector{Int32}}}(undef,Dc+1)
  d_to_patch_to_dpfaces[Dc+1] = patch_cells
  PatchTopology(topo,d_to_patch_to_dpfaces)
end

function PatchTopology(topo::GridTopology{Dc},Df=Dc) where Dc
  patch_cells = get_faces(topo,Df,Dc)
  PatchTopology(topo,patch_cells)
end

function PatchTopology(model::DiscreteModel,args...)
  topo = get_grid_topology(model)
  PatchTopology(topo,args...)
end

function num_patches(ptopo::PatchTopology)
  length(get_patch_cells(ptopo))
end

num_cells(ptopo::PatchTopology{Dc}) where Dc = num_faces(ptopo,Dc)
num_faces(ptopo::PatchTopology,d) = length(get_patch_faces(ptopo,d).data)

get_patch_cells(ptopo::PatchTopology{Dc}) where Dc = get_patch_faces(ptopo,Dc)
get_patch_facets(ptopo::PatchTopology{Dc}) where Dc = get_patch_faces(ptopo,Dc-1)

"""
    get_patch_faces(ptopo::PatchTopology,d::Integer) -> patch_faces
"""
function get_patch_faces(ptopo::PatchTopology,d::Integer)
  if isassigned(ptopo.d_to_patch_to_dfaces,d+1)
    patch_faces = ptopo.d_to_patch_to_dfaces[d+1]
  else
    patch_faces = generate_patch_faces(ptopo,d)
    ptopo.d_to_patch_to_dfaces[d+1] = patch_faces
  end
  return patch_faces
end

"""
    patch_reindex(ptopo::PatchTopology{Dc},face_to_data,Df=Dc) -> pface_to_data
"""
function patch_reindex(ptopo::PatchTopology{Dc},face_to_data,Df=Dc) where Dc
  patch_faces = get_patch_faces(ptopo,Df)
  pface_to_data = lazy_map(Reindex(face_to_data),patch_faces.data)
  return pface_to_data
end

"""
    patch_extend(PD::PatchTopology{Dc},patch_to_data,Df=Dc) -> pface_to_data
"""
function patch_extend(ptopo::PatchTopology{Dc},patch_to_data,Df=Dc) where Dc
  pface_to_patch = fill(0,num_faces(ptopo,Df))
  patch_offsets = get_patch_faces(ptopo,Df).ptrs
  for patch in 1:num_patches(ptopo)
    pface_to_patch[patch_offsets[patch]:patch_offsets[patch+1]-1] .= patch
  end

  pface_to_data = lazy_map(Reindex(patch_to_data),pface_to_patch)
  return pface_to_data
end

function generate_patch_faces(ptopo::PatchTopology{Dc},Df) where Dc
  cell_to_faces = get_faces(ptopo.topo,Dc,Df)

  n_patches = num_patches(ptopo)
  patch_cells = get_patch_cells(ptopo)
  
  pfaces = SortedSet{Int32}()
  ptrs = zeros(Int32,n_patches+1)
  for patch in 1:n_patches
    cells = view(patch_cells,patch)
    for cell in cells
      push!(pfaces, view(cell_to_faces,cell)...)
    end
    ptrs[patch+1] += length(pfaces)
    empty!(pfaces)
  end
  Arrays.length_to_ptrs!(ptrs)

  data = zeros(Int,ptrs[end]-1)
  empty!(pfaces)
  for patch in 1:n_patches
    cells = view(patch_cells,patch)
    for cell in cells
      push!(pfaces, view(cell_to_faces,cell)...)
    end
    data[ptrs[patch]:ptrs[patch+1]-1] = collect(pfaces)
    empty!(pfaces)
  end

  patch_faces = Table(data,ptrs)
  return patch_faces
end

function get_patch_boundary_info(ptopo::PatchTopology{Dc}) where Dc
  Df = Dc - 1 # TODO: Make general
  topo = ptopo.topo
  patch_cells = get_patch_cells(ptopo)
  patch_faces = get_patch_faces(ptopo,Df)
  face_to_cells = get_faces(topo,Df,Dc)

  n_pfaces = num_faces(ptopo,Dc-1)
  pface_to_isboundary = fill(false,n_pfaces)
  pface_to_lcell = fill(Int8(1),n_pfaces)

  k = 1
  for patch in 1:num_patches(ptopo)
    pcells = view(patch_cells,patch)
    pfaces = view(patch_faces,patch)
    for pface in pfaces
      pface_cells = view(face_to_cells,pface)
      lcells = findall(c -> c in pcells, pface_cells)
      if isone(length(lcells))
        pface_to_isboundary[k] = true
        pface_to_lcell[k] = Int8(lcells[1])
        k += 1
      end
    end
  end

  return pface_to_isboundary, pface_to_lcell
end

function get_pface_to_patch(ptopo::PatchTopology,Df)
  patch_faces = get_patch_faces(ptopo,Df)
  patch_pfaces = Table(Base.OneTo(num_faces(ptopo,Df)),patch_faces.ptrs)
  return Arrays.flatten_partition(patch_pfaces)
end

# PatchTriangulation

"""
    struct PatchTriangulation{Dc,Dp} <: Triangulation{Dc,Dp}
      ...
    end

Wrapper around a Triangulation, for patch-based assembly.
"""
struct PatchTriangulation{Dc,Dp,A,B,C} <: Triangulation{Dc,Dp}
  trian :: A
  ptopo :: B
  tface_to_pface :: C
  function PatchTriangulation(
    trian::Triangulation{Dc,Dp},
    ptopo::PatchTopology,
    tface_to_pface::AbstractVector{<:Integer}
  ) where {Dc,Dp}
    A, B, C = typeof(trian), typeof(ptopo), typeof(tface_to_pface)
    new{Dc,Dp,A,B,C}(trian,ptopo,tface_to_pface)
  end
end

# Triangulation API

get_background_model(t::PatchTriangulation) = get_background_model(t.trian)
get_grid(t::PatchTriangulation) = get_grid(t.trian)
get_glue(t::PatchTriangulation,::Val{D}) where D = get_glue(t.trian,Val(D))
get_facet_normal(trian::PatchTriangulation) = get_facet_normal(trian.trian)

# Constructors

function PatchTriangulation(model::DiscreteModel{Dc},ptopo::PatchTopology{Dc}) where Dc
  patch_cells = get_patch_cells(ptopo)
  trian = Triangulation(model,patch_cells.data)
  tcell_to_pcell = Base.OneTo(num_faces(ptopo,Dc))
  return PatchTriangulation(trian,ptopo,tcell_to_pcell)
end

function PatchBoundaryTriangulation(model::DiscreteModel{Dc},ptopo::PatchTopology{Dc};tags=nothing) where Dc
  patch_faces = get_patch_faces(ptopo,Dc-1)
  pface_to_face = patch_faces.data

  # Create p-face mask
  if !isnothing(tags)
    labels = get_face_labeling(model)
    face_to_mask = get_face_mask(labels,tags,Dc-1)
    pface_to_mask = lazy_map(Reindex(face_to_mask),pface_to_face)
  else
    pface_to_mask = Fill(true,num_faces(ptopo,Dc-1))
  end

  # Merge with boundary mask
  pface_to_isboundary, pface_to_lcell = get_patch_boundary_info(ptopo)
  tface_to_pface = findall(lazy_map(&,pface_to_isboundary,pface_to_mask))
  tface_to_lcell = pface_to_lcell[tface_to_pface]
  tface_to_face  = pface_to_face[tface_to_pface]

  # Create overlapping boundary triangulation
  topo = get_grid_topology(model)
  cell_grid = get_grid(model)
  face_grid = view(Grid(ReferenceFE{Dc-1},model),tface_to_face)
  face_glue = OverlappingFaceToCellGlue(topo,cell_grid,face_grid,tface_to_face,tface_to_lcell)
  face_trian = BodyFittedTriangulation(model,face_grid,tface_to_face)
  trian = BoundaryTriangulation(face_trian,face_glue)

  return PatchTriangulation(trian,ptopo,tface_to_pface)
end
