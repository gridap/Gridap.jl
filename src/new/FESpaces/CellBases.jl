# Abstract types

"""
"""
abstract type CellBasis <: CellField end

"""
"""
function TrialStyle(::Type{<:CellBasis})
  Val{false}()
end

TrialStyle(cb::CellBasis) = TrialStyle(typeof(cb))

"""
"""
is_trial(cb::CellBasis) = is_trial(typeof(cb))

is_test(cb::CellBasis) = ! is_trial(cb)

function is_trial(::Type{T}) where T<:CellBasis
  get_val_parameter(TrialStyle(T))
end

function is_test(::Type{T}) where T<:CellBasis
  ! is_trial(T)
end

"""
"""
function test_cell_basis(cb::CellBasis,args...; kwargs...)
  test_cell_field(cb,args...;kwargs...)
end

"""
"""
abstract type CellMatrixField <: CellField end

"""
"""
function test_cell_matrix_field(cb::CellMatrixField,args...; kwargs...)
  test_cell_field(cb,args...;kwargs...)
end

# Define how the metadata is preserved

function similar_cell_field(cf::CellBasis,a::AbstractArray)
  cm = get_cell_map(cf)
  trial_style = TrialStyle(cf)
  GenericCellBasis(trial_style,a,cm)
end

function similar_cell_field(cf::CellMatrixField,a::AbstractArray)
  cm = get_cell_map(cf)
  GenericCellMatrixField(a,cm)
end

function similar_cell_field(a::CellBasis,b::CellField,v::AbstractArray)
  similar_cell_field(a,v)
end

function similar_cell_field(a::CellField,b::CellBasis,v::AbstractArray)
  similar_cell_field(b,v)
end

function similar_cell_field(::CellBasis,::CellMatrixField,::AbstractArray)
  @notimplemented
end

function similar_cell_field(::CellMatrixField,::CellBasis,a::AbstractArray)
  @notimplemented
end

function similar_cell_field(a::CellMatrixField,b::CellField,v::AbstractArray)
  similar_cell_field(a,v)
end

function similar_cell_field(a::CellField,b::CellMatrixField,v::AbstractArray)
  similar_cell_field(b,v)
end

function similar_cell_field(a::CellMatrixField,b::CellMatrixField,v::AbstractArray)
  similar_cell_field(a,v)
end

function similar_cell_field(a::CellBasis,b::CellBasis,v::AbstractArray)
  _similar_cell_basis(v,a,b,TrialStyle(a),TrialStyle(b))
end

function _similar_cell_basis(v,a,b,a_trial::Val{T},b_trial::Val{T}) where T
  @notimplemented
end

function _similar_cell_basis(v,a,b,a_trial::Val{false},b_trial::Val{true})
  _similar_cell_basis_test_trial(v,a,b)
end

function _similar_cell_basis(v,a,b,a_trial::Val{true},b_trial::Val{false})
  _similar_cell_basis_test_trial(v,b,a)
end

function _similar_cell_basis_test_trial(v,a,b)
  cm = get_cell_map(a)
  GenericCellMatrixField(v,cm)
end

# Define operations

function operate_cell_field(op,a::CellBasis,b::CellMatrixField)
  @notimplemented
end

function operate_cell_field(op,a::CellMatrixField,b::CellBasis)
  @notimplemented
end

function operate_cell_field(op,a::CellBasis,b::CellBasis)
  _operate_cell_basis(op,a,b,TrialStyle(a),TrialStyle(b))
end

function _operate_cell_basis(op,a,b,atrial::Val{T},btrial::Val{T}) where T
  s = "Not implemented feature: it is unlikely that we need to operate between"
  s *= " two test (resp. two trial) CellBasis objects"
  @notimplemented s
end

function _operate_cell_basis(op,a,b,atrial::Val{false},btrial::Val{true})
  _operate_cell_basis_test_trial(op,a,b)
end

function _operate_cell_basis(op,a,b,atrial::Val{true},btrial::Val{false})
  _operate_cell_basis_test_trial(op,b,a)
end

function _operate_cell_basis_test_trial(op,a,b)
  operate_cell_field_default(op,a,b)
end

function operate_cell_field(op,a::CellMatrixField,b::CellMatrixField)
  _operate_cell_matrix_field(op,a,b)
end

function operate_cell_field(op,a::CellMatrixField,b::CellField)
  _operate_cell_matrix_field(op,a,b)
end

function operate_cell_field(op,a::CellField,b::CellMatrixField)
  _operate_cell_matrix_field(op,a,b)
end

function _operate_cell_matrix_field(op,a,b)
  aa = get_array(a)
  ab = get_array(b)
  k = bcast(op)
  c = apply_to_field_array(UnimplementedField,k,aa,ab)
  similar_cell_field(a,b,c)
end

# Concrete types

"""
"""
struct GenericCellBasis{T,A,B} <: CellBasis
  trial_style::Val{T}
  array::A
  cell_map::B
end

"""
"""
function GenericCellBasis(array::AbstractArray,cell_map::AbstractArray)
  trial_style = Val{false}()
  GenericCellBasis(trial_style,array,cell_map)
end

get_array(a::GenericCellBasis) = a.array

get_cell_map(a::GenericCellBasis) = a.cell_map

function TrialStyle(::Type{<:GenericCellBasis{T}}) where T
  Val{T}()
end

"""
"""
struct GenericCellMatrixField{A,B} <: CellMatrixField
  array::A
  cell_map::B
end

get_array(a::GenericCellMatrixField) = a.array

get_cell_map(a::GenericCellMatrixField) = a.cell_map

