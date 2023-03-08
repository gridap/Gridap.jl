
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
                          n2o_cell_to_child_id::AbstractVector{<:Integer},
                          refinement_rules::AbstractVector{<:RefinementRule})
    Dc = length(n2o_faces_map)-1
    is_refined    = select_refined_cells(n2o_faces_map[Dc+1])
    o2n_faces_map = get_o2n_faces_map(n2o_faces_map[Dc+1])

    GT = all(is_refined) ? RefinementGlue : MixedGlue

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
  @notimplemented
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

function get_new_cell_refinement_rules(g::AdaptivityGlue)
  old_rrules = g.refinement_rules
  n2o_faces_map = g.n2o_faces_map[end]
  return lazy_map(Reindex(old_rrules),n2o_faces_map)
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
