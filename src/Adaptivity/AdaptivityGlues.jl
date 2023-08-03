
abstract type AdaptivityGlueType end
struct RefinementGlue <: AdaptivityGlueType end
struct MixedGlue <: AdaptivityGlueType end

"""
Adaptivity glue between two nested triangulations:

- `n2o_faces_map`       : Given a new face gid, returns 
        A) if fine, the gid of the old face containing it.
        B) if coarse, the gids of its children (in child order)
- `n2o_cell_to_child_id`: Given a new cell gid, returns 
        A) if fine, the local child id within the (old) coarse cell containing it.
        B) if coarse, -1
- `refinement_rules`    : RefinementRule used for each coarse cell.
"""
struct AdaptivityGlue{GT,Dc,A,B,C,D,E} <: GridapType
  n2o_faces_map        :: A
  n2o_cell_to_child_id :: B
  refinement_rules     :: C
  o2n_faces_map        :: D
  is_refined           :: E

  function AdaptivityGlue(n2o_faces_map::Vector{<:Union{AbstractVector{<:Integer},Table{<:Integer}}},
                          n2o_cell_to_child_id::Union{AbstractVector{<:Integer},Table{<:Integer}},
                          refinement_rules::AbstractVector{<:RefinementRule})
    Dc = length(n2o_faces_map)-1
    is_refined    = select_refined_cells(n2o_faces_map[Dc+1])
    o2n_faces_map = get_o2n_faces_map(n2o_faces_map[Dc+1])

    GT = all(is_refined) ? RefinementGlue : MixedGlue
    if isa(GT(),RefinementGlue)
      @assert isa(n2o_faces_map,Vector{<:AbstractVector{<:Integer}})
      @assert isa(n2o_cell_to_child_id,AbstractVector{<:Integer})
    else
      @assert isa(n2o_faces_map,Vector{<:Table{<:Integer}})
      @assert isa(n2o_cell_to_child_id,Table{<:Integer})
    end

    A = typeof(n2o_faces_map)
    B = typeof(n2o_cell_to_child_id)
    C = typeof(refinement_rules)
    D = typeof(o2n_faces_map)
    E = typeof(is_refined)
    new{GT,Dc,A,B,C,D,E}(n2o_faces_map,n2o_cell_to_child_id,refinement_rules,o2n_faces_map,is_refined)
  end
end

function AdaptivityGlue(n2o_faces_map,n2o_cell_to_child_id,refinement_rule::RefinementRule)
  n_old  = maximum(n2o_faces_map[end])
  rrules = Fill(refinement_rule,n_old)
  return AdaptivityGlue(n2o_faces_map,n2o_cell_to_child_id,rrules)
end

select_refined_cells(n2o_cell_map::Vector) = Fill(true,length(n2o_cell_map))
select_refined_cells(n2o_cell_map::Table) = map(x -> length(x) == 1, n2o_cell_map)

"""
For each fine cell, returns Φ st. x_coarse = ϕ(x_fine)
"""
function get_n2o_reference_coordinate_map(g::AdaptivityGlue{RefinementGlue})
  rrules    = get_new_cell_refinement_rules(g)
  cell_maps = lazy_map(Geometry.get_cell_map,rrules)
  return lazy_map(getindex,cell_maps,g.n2o_cell_to_child_id)
end

function get_n2o_reference_coordinate_map(g::AdaptivityGlue{MixedGlue})
  @notimplemented
end

"""
  Given a map from new cells to old cells, computes the inverse map.
  In the general case (refinement+coarsening), the n2o map is a Table, 
  but the algorithm is optimized for Vectors (refinement only).
"""
function get_o2n_faces_map(ncell_to_ocell::Table{T}) where {T<:Integer}
  nC = maximum(ncell_to_ocell.data)
  
  ptrs = fill(0,nC+1)
  for ccell in ncell_to_ocell.data
    ptrs[ccell+1] += 1
  end
  Arrays.length_to_ptrs!(ptrs)

  data = Vector{Int}(undef,ptrs[end]-1)
  for fcell = 1:length(ncell_to_ocell.ptrs)-1
    for j = ncell_to_ocell.ptrs[fcell]:ncell_to_ocell.ptrs[fcell+1]-1
      ccell = ncell_to_ocell.data[j]
      data[ptrs[ccell]] = fcell
      ptrs[ccell] += 1
    end
  end
  Arrays.rewind_ptrs!(ptrs)

  ocell_to_ncell = Table(data,ptrs)
  return ocell_to_ncell
end

function get_o2n_faces_map(ncell_to_ocell::Vector{T}) where {T<:Integer}
  (length(ncell_to_ocell) == 0) && (return Table(T[],T[]))

  nC = maximum(ncell_to_ocell)
  nF = length(ncell_to_ocell)

  ptrs = fill(0,nC+1)
  for iF in 1:nF
    iC = ncell_to_ocell[iF]
    ptrs[iC+1] += 1
  end
  length_to_ptrs!(ptrs)

  cnts = fill(0,nC)
  data = fill(zero(T),ptrs[end])
  for iF in 1:nF
    iC = ncell_to_ocell[iF]
    data[ptrs[iC]+cnts[iC]] = iF
    cnts[iC] += 1
  end

  ocell_to_ncell = Table(data,ptrs)
  return ocell_to_ncell
end

function get_new_cell_refinement_rules(g::AdaptivityGlue{<:RefinementGlue})
  old_rrules = g.refinement_rules
  n2o_faces_map = g.n2o_faces_map[end]
  return lazy_map(Reindex(old_rrules),n2o_faces_map)
end

function get_new_cell_refinement_rules(g::AdaptivityGlue{<:MixedGlue})
  old_rrules = g.refinement_rules
  n2o_faces_map = g.n2o_faces_map[end]
  new_idx = lazy_map(Reindex(n2o_faces_map.data),n2o_faces_map.ptrs[1:end-1])
  return lazy_map(Reindex(old_rrules), new_idx)
end

function get_old_cell_refinement_rules(g::AdaptivityGlue)
  return g.refinement_rules
end

# Data re-indexing

function f2c_reindex(fine_data,g::AdaptivityGlue)
  ccell_to_fcell = g.o2n_faces_map
  return _reindex(fine_data,ccell_to_fcell)
end

function c2f_reindex(coarse_data,g::AdaptivityGlue{GT,Dc}) where {GT,Dc}
  fcell_to_ccell = g.n2o_faces_map[Dc+1]
  return _reindex(coarse_data,fcell_to_ccell)
end

f2c_reindex(a::CellDatum,g::AdaptivityGlue) = f2c_reindex(CellData.get_data(a),g::AdaptivityGlue)
c2f_reindex(a::CellDatum,g::AdaptivityGlue) = c2f_reindex(CellData.get_data(a),g::AdaptivityGlue)

function _reindex(data,idx::Table)
  m = Reindex(data)
  return Table(lazy_map(m,idx.data),idx.ptrs)
end

function _reindex(data,idx::Vector)
  m = Reindex(data)
  return lazy_map(m,idx)
end


# New to old face glues 

function get_d_to_fface_to_cface(::AdaptivityGlue,::GridTopology,::GridTopology)
  @notimplemented
end

"
 For each child/fine face, returns the parent/coarse face containing it. The parent 
 face might have higher dimension. 

 Returns two arrays: 
  - [dimension][fine face gid] -> coarse parent face gid
  - [dimension][fine face gid] -> coarse parent face dimension
"
function get_d_to_fface_to_cface(glue::AdaptivityGlue{<:RefinementGlue},
                                 ctopo::GridTopology{Dc},
                                 ftopo::GridTopology{Dc}) where Dc

  # Local data for each coarse cell, at the RefinementRule level
  rrules = Adaptivity.get_old_cell_refinement_rules(glue)
  ccell_to_d_to_faces = lazy_map(rr->map(d->Geometry.get_faces(get_grid_topology(rr.ref_grid),Dc,d),0:Dc),rrules)
  ccell_to_d_to_fface_to_parent_face = lazy_map(get_d_to_face_to_parent_face,rrules)

  # Global data, concerning the complete meshes
  ccell_to_fcell = glue.o2n_faces_map
  d_to_ccell_to_cface = map(d->Geometry.get_faces(ctopo,Dc,d),0:Dc)
  d_to_fcell_to_fface = map(d->Geometry.get_faces(ftopo,Dc,d),0:Dc)

  d_to_fface_to_cface = [fill(Int32(0),num_faces(ftopo,d)) for d in 0:Dc]
  d_to_fface_to_cface_dim = [fill(Int32(0),num_faces(ftopo,d)) for d in 0:Dc]

  # For each coarse cell
  for ccell in 1:num_cells(ctopo)
    local_d_to_fface_to_parent_face,
      local_d_to_fface_to_parent_dim = ccell_to_d_to_fface_to_parent_face[ccell]
    # For each fine subcell:
    # child_id -> Local Id of the fine cell within the refinement rule (ccell)
    for (child_id,fcell) in enumerate(ccell_to_fcell[ccell])
      # For each fine face on the fine subcell: 
      # d     -> Dimension of the fine face
      # iF    -> Local Id of the fine face within the fine cell
      # fface -> Global Id of the fine face 
      for d in 0:Dc
        for (iF,fface) in enumerate(d_to_fcell_to_fface[d+1][fcell])
          # Local Id of the fine face within the refinement rule
          fface_child_id = ccell_to_d_to_faces[ccell][d+1][child_id][iF]
          # Local Id of the coarse parent face within the coarse cell
          parent    = local_d_to_fface_to_parent_face[d+1][fface_child_id]

          # Global Id of the coarse parent face, and it's dimension
          cface_dim = local_d_to_fface_to_parent_dim[d+1][fface_child_id]
          cface     = d_to_ccell_to_cface[cface_dim+1][ccell][parent]
          d_to_fface_to_cface[d+1][fface] = cface
          d_to_fface_to_cface_dim[d+1][fface] = cface_dim
        end
      end
    end
  end

  return (d_to_fface_to_cface,d_to_fface_to_cface_dim)
end

# FaceLabeling refinement

function _refine_face_labeling(coarse_labeling::FaceLabeling,
                               glue  :: AdaptivityGlue,
                               ctopo :: GridTopology,
                               ftopo :: GridTopology)
  d_to_fface_to_cface,
    d_to_fface_to_cface_dim = get_d_to_fface_to_cface(glue,ctopo,ftopo)

  return _refine_face_labeling(coarse_labeling,d_to_fface_to_cface,d_to_fface_to_cface_dim)
end

function _refine_face_labeling(coarse_labeling::FaceLabeling,
                               d_to_fface_to_cface,
                               d_to_fface_to_cface_dim)
  tag_to_name = copy(coarse_labeling.tag_to_name)
  tag_to_entities = copy(coarse_labeling.tag_to_entities)
  
  Dc = num_dims(coarse_labeling)
  d_to_dface_to_entity = Vector{Vector{Int32}}(undef,Dc+1)
  for d in 0:Dc
    nF = length(d_to_fface_to_cface[d+1])
    dface_to_entity = Vector{Int32}(undef,nF)
  
    for fface in 1:nF
      cface = d_to_fface_to_cface[d+1][fface]
      cface_dim = d_to_fface_to_cface_dim[d+1][fface]

      cface_entity = coarse_labeling.d_to_dface_to_entity[cface_dim+1][cface]
      dface_to_entity[fface] = cface_entity
    end
  
    d_to_dface_to_entity[d+1] = dface_to_entity
  end
  
  return Geometry.FaceLabeling(d_to_dface_to_entity,tag_to_entities,tag_to_name)  
end