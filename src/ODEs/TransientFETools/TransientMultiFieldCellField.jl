################################
# TransientMultiFieldCellField #
################################
struct TransientMultiFieldCellField{A} <: TransientCellField
  cellfield::A
  derivatives::Tuple
  transient_single_fields::Vector{<:TransientCellField} # used to iterate
end

const MultiFieldTypes = Union{MultiFieldCellField,MultiFieldFEFunction}

function TransientCellField(multi_field::MultiFieldTypes, derivatives::Tuple)
  _flat = _to_transient_single_fields(multi_field, derivatives)
  TransientMultiFieldCellField(multi_field, derivatives, _flat)
end

# CellField interface
function get_data(f::TransientMultiFieldCellField)
  s = """
  Function get_data is not implemented for TransientMultiFieldCellField at this moment.
  You need to extract the individual fields and then evaluate them separately.

  If ever implement this, evaluating a `MultiFieldCellField` directly would provide,
  at each evaluation point, a tuple with the value of the different fields.
  """
  @notimplemented s
end

get_triangulation(f::TransientMultiFieldCellField) = get_triangulation(f.cellfield)

DomainStyle(::Type{TransientMultiFieldCellField{A}}) where {A} = DomainStyle(A)

gradient(f::TransientMultiFieldCellField) = gradient(f.cellfield)

∇∇(f::TransientMultiFieldCellField) = ∇∇(f.cellfield)

function change_domain(f::TransientMultiFieldCellField, trian::Triangulation,
  target_domain::DomainStyle)
  change_domain(f.cellfield, trian, target_domain)
end

# MultiField interface
num_fields(f::TransientMultiFieldCellField) = length(f.cellfield)

# Base interface
function getindex(f::TransientMultiFieldCellField, index::Integer)
  sub_cellfield = f.cellfield[index]

  sub_derivatives = ()
  for derivative in f.derivatives
    sub_derivative = derivative[index]
    sub_derivatives = (sub_derivatives..., sub_derivative)
  end

  TransientSingleFieldCellField(sub_cellfield, sub_derivatives)
end

function getindex(f::TransientMultiFieldCellField,
  indices::AbstractVector{<:Integer})
  sub_cellfield = MultiFieldCellField(
    f.cellfield[indices],
    DomainStyle(f.cellfield)
  )

  sub_derivatives = ()
  for derivative in f.derivatives
    sub_derivative = MultiFieldCellField(
      derivative[indices],
      DomainStyle(derivative)
    )
    sub_derivatives = (sub_derivatives..., sub_derivative)
  end

  _sub_flat = _to_transient_single_fields(sub_cellfield, sub_derivatives)
  TransientMultiFieldCellField(sub_cellfield, sub_derivatives, _sub_flat)
end

function iterate(f::TransientMultiFieldCellField)
  iterate(f.transient_single_fields)
end

function iterate(f::TransientMultiFieldCellField, state)
  iterate(f.transient_single_fields, state)
end

# Time derivatives
function ∂t(f::TransientMultiFieldCellField)
  cellfield, derivatives = first_and_tail(f.derivatives)

  transient_single_field_derivatives = TransientCellField[]
  for transient_single_field in f.transient_single_fields
    push!(transient_single_field_derivatives, ∂t(transient_single_field))
  end

  TransientMultiFieldCellField(
    cellfield, derivatives,
    transient_single_field_derivatives
  )
end

∂tt(f::TransientMultiFieldCellField) = ∂t(∂t(f))

#########
# Utils #
#########
function _to_transient_single_fields(multi_field, derivatives)
  transient_single_fields = TransientCellField[]

  for index in 1:num_fields(multi_field)
    single_field = multi_field[index]

    single_derivatives = ()
    for derivative in derivatives
      single_derivatives = (single_derivatives..., derivative[index])
    end

    transient_single_field = TransientSingleFieldCellField(
      single_field,
      single_derivatives
    )
    push!(transient_single_fields, transient_single_field)
  end

  transient_single_fields
end
