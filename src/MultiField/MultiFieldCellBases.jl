
struct CellBasisWithFieldID{S,R} <: CellBasis
  trial_style::Val{S}
  cell_basis::CellBasis
  field_id::Int
  ref_style::Val{R}
  function CellBasisWithFieldID(cell_basis::CellBasis,field_id::Integer)
    trial_style = TrialStyle(cell_basis)
    S = get_val_parameter(trial_style)
    ref_style = RefStyle(cell_basis)
    R = get_val_parameter(ref_style)
    new{S,R}(trial_style,cell_basis,field_id,ref_style)
  end
end

get_array(cb::CellBasisWithFieldID) = get_array(cb.cell_basis)

get_cell_map(cb::CellBasisWithFieldID) = get_cell_map(cb.cell_basis)

TrialStyle(::Type{CellBasisWithFieldID{S,R}}) where {S,R} = Val{S}()

RefStyle(::Type{CellBasisWithFieldID{S,R}}) where {S,R} = Val{R}()

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
  @assert is_in_ref_space(a) == is_in_ref_space(b)
  CellMatrixFieldWithFieldIds(cell_matrix_field,field_id_rows,field_id_cols,RefStyle(a))
end

function operate(op,a::CellBasisWithFieldID,b::CellBasisWithFieldID)
  _operate_cell_basis_with_field_id(op,a,b,TrialStyle(a),TrialStyle(b))
end

function _operate_cell_basis_with_field_id(op,a,b,atrial::Val{T},btrial::Val{T}) where T
  if a.field_id == b.field_id
    _operate_cell_basis(op,a,b,TrialStyle(a),TrialStyle(b))
  else
    _a = BlockTracker(a)
    _b = BlockTracker(b)
    op(_a,_b)
  end
end

function BlockTracker(a::CellBasisWithFieldID)
  blocks = (a,)
  block_ids = [(a.field_id,)]
  BlockTracker(blocks,block_ids)
end

operate(op,a::BlockTracker,b::CellBasisWithFieldID) = operate(op,a,BlockTracker(b))

operate(op,a::CellBasisWithFieldID,b::BlockTracker) = operate(op,BlockTracker(a),b)

function _operate_cell_basis_with_field_id(op,a,b,atrial,btrial)
  _operate_cell_basis(op,a,b,TrialStyle(a),TrialStyle(b))
end

struct CellMatrixFieldWithFieldIds{R} <: CellMatrixField
  cell_matrix_field::CellMatrixField
  field_id_rows::Int
  field_id_cols::Int
  ref_style::Val{R}
end

RefStyle(::Type{CellMatrixFieldWithFieldIds{R}}) where R = Val{R}()

get_array(a::CellMatrixFieldWithFieldIds) = get_array(a.cell_matrix_field)

get_cell_map(a::CellMatrixFieldWithFieldIds) = get_cell_map(a.cell_matrix_field)

function _have_same_field_ids(a::CellMatrixFieldWithFieldIds,b::CellMatrixFieldWithFieldIds)
  (a.field_id_rows == b.field_id_rows) && (a.field_id_cols == b.field_id_cols)
end

function similar_object(cf::CellMatrixFieldWithFieldIds,a::AbstractArray)
  cell_matrix_field = similar_object(cf.cell_matrix_field,a)
  field_id_rows = cf.field_id_rows
  field_id_cols = cf.field_id_cols
  CellMatrixFieldWithFieldIds(cell_matrix_field,field_id_rows,field_id_cols,RefStyle(cf))
end

function similar_object(a::CellMatrixFieldWithFieldIds,b::CellMatrixFieldWithFieldIds,v::AbstractArray)
  @assert _have_same_field_ids(a,b)
  similar_object(a,v)
end

function operate(op,a::CellMatrixFieldWithFieldIds,b::CellMatrixFieldWithFieldIds)
  if _have_same_field_ids(a,b)
    _operate_cell_matrix_field(op,a,b)
  else
    _a = BlockTracker(a)
    _b = BlockTracker(b)
    op(_a,_b)
  end
end

function BlockTracker(a::CellMatrixFieldWithFieldIds)
  blocks = (a,)
  block_ids = [(a.field_id_rows,a.field_id_cols),]
  BlockTracker(blocks,block_ids)
end

operate(op,a::BlockTracker,b::CellMatrixFieldWithFieldIds) = operate(op,a,BlockTracker(b))

operate(op,a::CellMatrixFieldWithFieldIds,b::BlockTracker) = operate(op,BlockTracker(a),b)

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

# Evaluation

function evaluate(bs::AbstractArray{<:CellBasisWithFieldID},q::AbstractArray)
  f = b -> evaluate(b,q)
  blocks = tuple(map(f,bs)...)
  coordinates = [ (i,1) for i in 1:length(bs)]
  MultiFieldCellArray(blocks,coordinates)
end

# Integration

function integrate(cb::CellBasisWithFieldID,trian::Triangulation,quad::CellQuadrature)
  r = integrate(cb.cell_basis,trian,quad)
  bloks = (r,)
  block_ids = [(cb.field_id,),]
  MultiFieldCellArray(bloks,block_ids)
end

function integrate(cm::CellMatrixFieldWithFieldIds,trian::Triangulation,quad::CellQuadrature)
  r = integrate(cm.cell_matrix_field,trian,quad)
  bloks = (r,)
  block_ids = [(cm.field_id_rows, cm.field_id_cols),]
  MultiFieldCellArray(bloks,block_ids)
end

function integrate(cb::BlockTracker,trian::Triangulation,quad::CellQuadrature)
  f = (b) -> integrate(get_array(b),trian,quad)
  blocks = map(f,cb.blocks)
  block_ids = cb.block_ids
  MultiFieldCellArray(blocks,block_ids)
end

struct MultiCellBasis{S,R} <: GridapType
  blocks::Vector{CellBasisWithFieldID{S}}
  function MultiCellBasis(blocks::Vector{<:CellBasis})
    cb = first(blocks)
    S = is_trial(cb)
    R = is_in_ref_space(cb)
    new_blocks = CellBasisWithFieldID{S}[]
    for (field_id, cell_basis) in enumerate(blocks)
      @assert is_trial(cell_basis) == S "All the provided bases must be of the same type (trial or test)"
      @assert is_in_ref_space(cell_basis) == R "All the provided bases must be defined in the same space (reference or physical)"
      block = CellBasisWithFieldID(cell_basis,field_id)
      push!(new_blocks,block)
    end
    new{S,R}(new_blocks)
  end
end

FECellBasisStyle(::Type{<:MultiCellBasis}) = Val{true}()

TrialStyle(::Type{MultiCellBasis{S,R}}) where {S,R} = Val{S}()

RefStyle(::Type{MultiCellBasis{S,R}}) where {S,R} = Val{R}()

num_fields(m::MultiCellBasis) = length(m.blocks)

Base.iterate(m::MultiCellBasis) = iterate(m.blocks)

Base.iterate(m::MultiCellBasis,state) = iterate(m.blocks,state)

Base.getindex(m::MultiCellBasis,field_id::Integer) = m.blocks[field_id]

function restrict(a::MultiCellBasis,trian::Triangulation)
  f = (ai) -> restrict(ai,trian)
  blocks = map(f,a.blocks)
  blocks
end

# Dirichlet related

function kernel_cache(k::DirichletVecKernel,mat::MultiFieldArray,vals::MultiFieldArray)
  vec = mat*vals
  cvec = CachedMultiFieldArray(vec)
  (vec, cvec)
end

@inline function apply_kernel!(cache,k::DirichletVecKernel,mat::MultiFieldArray,vals::MultiFieldArray)
  vec, cvec = cache
  _resize_for_mul!(cvec,mat,vals)
  _move_cached_arrays!(vec,cvec)
  mul!(vec,mat,vals)
  for vk in vec.blocks
    @inbounds for i in eachindex(vk)
      vk[i] = -vk[i]
    end
  end
  vec
end
