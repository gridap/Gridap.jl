
struct CellBasisWithFieldID{S} <: CellBasis
  trial_style::Val{S}
  cell_basis::CellBasis
  field_id::Int
  function CellBasisWithFieldID(cell_basis::CellBasis,field_id::Integer)
    trial_style = TrialStyle(cell_basis)
    S = get_val_parameter(trial_style)
    new{S}(trial_style,cell_basis,field_id)
  end
end

get_array(cb::CellBasisWithFieldID) = get_array(cb.cell_basis)

get_cell_map(cb::CellBasisWithFieldID) = get_cell_map(cb.cell_basis)

TrialStyle(::Type{CellBasisWithFieldID{S}}) where S = Val{S}()

function similar_object(cb::CellBasisWithFieldID,array::AbstractArray)
  cell_basis = similar_object(cb.cell_basis,array)
  field_id = cb.field_id
  CellBasisWithFieldID(cell_basis,field_id)
end

function similar_object(a::CellBasisWithFieldID,b::CellBasisWithFieldID,v::AbstractArray)
  _similar_cell_basis_with_field_id(v,a,b,TrialStyle(a),TrialStyle(b))
end

function _similar_cell_basis_with_field_id(v,a,b,a_trial::Val{T},b_trial::Val{T}) where T
  cell_basis = similar_object(a.cell_basis,b.cell_basis,v)
  @assert a.field_id == b.field_id
  field_id = a.field_id
  CellBasisWithFieldID(cell_basis,field_id)
end

function _similar_cell_basis_with_field_id(v,a,b,a_trial::Val{false},b_trial::Val{true})
  _similar_cell_basis_with_field_id_test_trial(v,a,b)
end

function _similar_cell_basis_with_field_id(v,a,b,a_trial::Val{true},b_trial::Val{false})
  _similar_cell_basis_with_field_id_test_trial(v,b,a)
end

function  _similar_cell_basis_with_field_id_test_trial(v,a,b)
  field_id_rows = a.field_id
  field_id_cols = b.field_id
  cell_matrix_field = similar_object(a.cell_basis,b.cell_basis,v)
  CellMatrixFieldWithFieldIds(cell_matrix_field,field_id_rows,field_id_cols)
end

function operate(op,a::CellBasisWithFieldID,b::CellBasisWithFieldID)
  _operate_cell_basis_with_field_id(op,a,b,TrialStyle(a),TrialStyle(b))
end

function _operate_cell_basis_with_field_id(op,a,b,atrial::Val{T},btrial::Val{T}) where T
  if a.field_id == b.field_id
    _operate_cell_basis(op,a,b,TrialStyle(a),TrialStyle(b))
  else
    _a = BlockTracker((a,),[(a.field_id,)])
    _b = BlockTracker((b,),[(b.field_id,)])
    op(_a,_b)
  end
end

function _operate_cell_basis_with_field_id(op,a,b,atrial,btrial)
  _operate_cell_basis(op,a,b,TrialStyle(a),TrialStyle(b))
end

struct CellMatrixFieldWithFieldIds <: CellMatrixField
  cell_matrix_field::CellMatrixField
  field_id_rows::Int
  field_id_cols::Int
end

get_array(a::CellMatrixFieldWithFieldIds) = get_array(a.cell_matrix_field)

get_cell_map(a::CellMatrixFieldWithFieldIds) = get_cell_map(a.cell_matrix_field)

function _have_same_field_ids(a::CellMatrixFieldWithFieldIds,b::CellMatrixFieldWithFieldIds)
  (a.field_id_rows == b.field_id_rows) && (a.field_id_cols == b.field_id_cols)
end

function similar_object(cf::CellMatrixFieldWithFieldIds,a::AbstractArray)
  cell_matrix_field = similar_object(cf.cell_matrix_field,a)
  field_id_rows = cf.field_id_rows
  field_id_cols = cf.field_id_cols
  CellMatrixFieldWithFieldIds(cell_matrix_field,field_id_rows,field_id_cols)
end

function similar_object(a::CellMatrixFieldWithFieldIds,b::CellMatrixFieldWithFieldIds,v::AbstractArray)
  @assert _have_same_field_ids(a,b)
  similar_object(a,v)
end

function operate(op,a::CellMatrixFieldWithFieldIds,b::CellMatrixFieldWithFieldIds)
  if _have_same_field_ids(a,b)
    _operate_cell_matrix_field(op,a,b)
  else
    _a = BlockTracker((a,),[(a.field_id_rows,a.field_id_cols)])
    _b = BlockTracker((b,),[(b.field_id_rows,b.field_id_cols)])
    op(_a,_b)
  end
end

# Restrictions

function restrict(cf::CellBasisWithFieldID,trian::Triangulation)
  r = restrict(cf.cell_basis,trian)
  _attach_field_id(r,cf.field_id)
end

function _attach_field_id(r::CellBasis,field_id)
  CellBasisWithFieldID(r,field_id)
end

function _attach_field_id(r::SkeletonCellBasis,field_id)
  l = CellBasisWithFieldID(r.left,field_id)
  r = CellBasisWithFieldID(r.right,field_id)
  SkeletonCellBasis(r.trial_style,l,r)
end

# Integration









#struct MultiCellBasis{S} <: GridapType
#  blocks::Vector{CellBasisWithFieldID{S}}
#  function MultiCellBasis(blocks::Vector{<:CellBasis})
#    cb = first(blocks)
#    S = is_trial(cb)
#    new_blocks = CellBasisWithFieldID{S}[]
#    for (field_id, cell_basis) in enumerate(blocks)
#      @assert is_trial(cell_basis) == S "All the provided bases need to be either test or trial"
#      block = CellBasisWithFieldID(cell_basis,field_id)
#      push!(new_blocks,block)
#    end
#    new{S}(new_blocks)
#  end
#end
#
#num_fields(m::MultiCellBasis) = length(m.blocks)
#
#Base.iterate(m::MultiCellBasis) = iterate(m.blocks)
#
#Base.iterate(m::MultiCellBasis,state) = iterate(m.blocks,state)
#
#Base.getindex(m::MultiCellBasis,field_id::Integer) = m.blocks[field_id]

#struct MultiCellArrayField{N} <: GridapType
#  blocks::Tuple
#  block_ids::Vector{NTuple{N,Int}}
#end
#
#function operate(op,a::MultiCellMatrixField)
#  f = (b)->operate(op,b)
#  new_blocks = map(f,a.blocks)
#  MultiCellMatrixField(new_blocks,a.block_ids)
#end
#
#function operate(op,a::MultiCellMatrixField,b)
#  f = (ai)->operate(op,ai,b)
#  new_blocks = map(f,a.blocks)
#  MultiCellMatrixField(new_blocks,a.block_ids)
#end
#
#function operate(op,a,b::MultiCellMatrixField)
#  f = (bi)->operate(op,a,bi)
#  new_blocks = map(f,b.blocks)
#  MultiCellMatrixField(new_blocks,b.block_ids)
#end
#
#function operate(op,a::MultiCellMatrixField,b::MultiCellMatrixField)
#  @assert op in (+,-)
#  a.blocks
#  f = (bi)->operate(op,bi)
#  b_blocks = map(f,b.blocks)
#  new_blocks = (a.blocks...,b_blocks...)
#  new_block_ids = vcat(a.block_ids,b.block_ids)
#  MultiCellMatrixField(new_blocks,new_block_ids)
#end


