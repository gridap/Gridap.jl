struct TransientMultiFieldCellField{A} <: TransientCellField
  cellfield::A
  derivatives::Tuple
  transient_single_fields::Vector{<:TransientCellField} # used to iterate
end

MultiFieldTypes = Union{MultiFieldCellField,MultiFieldFEFunction}

function TransientCellField(multi_field::MultiFieldTypes,derivatives::Tuple)
  transient_single_fields = _to_transient_single_fields(multi_field,derivatives)
  TransientMultiFieldCellField(multi_field,derivatives,transient_single_fields)
end

function get_data(f::TransientMultiFieldCellField)
  s = """
  Function get_data is not implemented for TransientMultiFieldCellField at this moment.
  You need to extract the individual fields and then evaluate them separatelly.

  If ever implement this, evaluating a `MultiFieldCellField` directly would provide,
  at each evaluation point, a tuple with the value of the different fields.
  """
  @notimplemented s
end

get_triangulation(f::TransientMultiFieldCellField) = get_triangulation(f.cellfield)
DomainStyle(::Type{TransientMultiFieldCellField{A}}) where A = DomainStyle(A)
num_fields(f::TransientMultiFieldCellField) = length(f.cellfield)
gradient(f::TransientMultiFieldCellField) = gradient(f.cellfield)
∇∇(f::TransientMultiFieldCellField) = ∇∇(f.cellfield)
change_domain(f::TransientMultiFieldCellField,trian::Triangulation,target_domain::DomainStyle) = change_domain(f.cellfield,trian,target_domain)

# Get single index
function Base.getindex(f::TransientMultiFieldCellField,ifield::Integer)
  single_field = f.cellfield[ifield]
  single_derivatives = ()
  for ifield_derivatives in f.derivatives
    single_derivatives = (single_derivatives...,getindex(ifield_derivatives,ifield))
  end
  TransientSingleFieldCellField(single_field,single_derivatives)
end

# Get multiple indices
function Base.getindex(f::TransientMultiFieldCellField,indices::Vector{<:Int})
  cellfield = MultiFieldCellField(f.cellfield[indices],DomainStyle(f.cellfield))
  derivatives = ()
  for derivative in f.derivatives
    derivatives = (derivatives...,MultiFieldCellField(derivative[indices],DomainStyle(derivative)))
  end
  transient_single_fields = _to_transient_single_fields(cellfield,derivatives)
  TransientMultiFieldCellField(cellfield,derivatives,transient_single_fields)
end

function _to_transient_single_fields(multi_field,derivatives)
  transient_single_fields = TransientCellField[]
  for ifield in 1:num_fields(multi_field)
    single_field = multi_field[ifield]
    single_derivatives = ()
    for ifield_derivatives in derivatives
      single_derivatives = (single_derivatives...,getindex(ifield_derivatives,ifield))
    end
    transient_single_field = TransientSingleFieldCellField(single_field,single_derivatives)
    push!(transient_single_fields,transient_single_field)
  end
  transient_single_fields
end

# Iterate functions
Base.iterate(f::TransientMultiFieldCellField)  = iterate(f.transient_single_fields)
Base.iterate(f::TransientMultiFieldCellField,state)  = iterate(f.transient_single_fields,state)

# Time derivative
function ∂t(f::TransientMultiFieldCellField)
  cellfield, derivatives = first_and_tail(f.derivatives)
  transient_single_field_derivatives = TransientCellField[]
  for transient_single_field in f.transient_single_fields
    push!(transient_single_field_derivatives,∂t(transient_single_field))
  end
  TransientMultiFieldCellField(cellfield,derivatives,transient_single_field_derivatives)
end

∂tt(f::TransientMultiFieldCellField) = ∂t(∂t(f))
