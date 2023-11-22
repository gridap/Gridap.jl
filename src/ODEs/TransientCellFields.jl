######################
# TransientCellField #
######################
"""
    abstract type TransientCellField <: CellField end

Transient CellField
"""
abstract type TransientCellField <: CellField end

# CellField interface
CellData.get_data(f::TransientCellField) = @abstractmethod

CellData.get_triangulation(f::TransientCellField) = @abstractmethod

CellData.DomainStyle(::Type{TransientCellField}) = @abstractmethod

function CellData.change_domain(
  f::TransientCellField, trian::Triangulation, target_domain::DomainStyle
)
  @abstractmethod
end

Fields.gradient(f::TransientCellField) = @abstractmethod

Fields.∇∇(f::TransientCellField) = @abstractmethod

#################################
# TransientSingleFieldCellField #
#################################
"""
    struct TransientSingleFieldCellField <: TransientCellField

Transient CellField for a SingleField FESpace
"""
struct TransientSingleFieldCellField{A} <: TransientCellField
  cellfield::A
  derivatives::Tuple # {Vararg{A,B} where B}
end

const SingleFieldTypes = Union{GenericCellField,SingleFieldFEFunction}

function TransientCellField(single_field::SingleFieldTypes, derivatives::Tuple)
  TransientSingleFieldCellField(single_field, derivatives)
end

# CellField interface
CellData.get_data(f::TransientSingleFieldCellField) = get_data(f.cellfield)

CellData.get_triangulation(f::TransientSingleFieldCellField) = get_triangulation(f.cellfield)

CellData.DomainStyle(::Type{<:TransientSingleFieldCellField{A}}) where {A} = DomainStyle(A)

function CellData.change_domain(f::TransientSingleFieldCellField, trian::Triangulation,
  target_domain::DomainStyle)
  change_domain(f.cellfield, trian, target_domain)
end

Fields.gradient(f::TransientSingleFieldCellField) = gradient(f.cellfield)

Fields.∇∇(f::TransientSingleFieldCellField) = ∇∇(f.cellfield)

# Skeleton-related operations
function Base.getproperty(f::TransientSingleFieldCellField, sym::Symbol)
  if sym in (:⁺, :plus, :⁻, :minus)
    derivatives = ()
    if sym in (:⁺, :plus)
      cellfield = CellFieldAt{:plus}(f.cellfield)
      for iderivative in f.derivatives
        derivatives = (derivatives..., CellFieldAt{:plus}(iderivative))
      end
    elseif sym in (:⁻, :minus)
      cellfield = CellFieldAt{:minus}(f.cellfield)
      for iderivative in f.derivatives
        derivatives = (derivatives..., CellFieldAt{:plus}(iderivative))
      end
    end
    return TransientSingleFieldCellField(cellfield, derivatives)
  else
    return getfield(f, sym)
  end
end

################################
# TransientMultiFieldCellField #
################################
"""
    struct TransientMultiFieldCellField <: TransientCellField

Transient CellField for a MultiField FESpace
"""
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
function CellData.get_data(f::TransientMultiFieldCellField)
  s = """
  Function get_data is not implemented for TransientMultiFieldCellField at this moment.
  You need to extract the individual fields and then evaluate them separately.

  If ever implement this, evaluating a `MultiFieldCellField` directly would provide,
  at each evaluation point, a tuple with the value of the different fields.
  """
  @notimplemented s
end

CellData.get_triangulation(f::TransientMultiFieldCellField) = get_triangulation(f.cellfield)

CellData.DomainStyle(::Type{TransientMultiFieldCellField{A}}) where {A} = DomainStyle(A)

function CellData.change_domain(f::TransientMultiFieldCellField, trian::Triangulation,
  target_domain::DomainStyle)
  change_domain(f.cellfield, trian, target_domain)
end

Fields.gradient(f::TransientMultiFieldCellField) = gradient(f.cellfield)

Fields.∇∇(f::TransientMultiFieldCellField) = ∇∇(f.cellfield)

# MultiField interface
MultiField.num_fields(f::TransientMultiFieldCellField) = length(f.cellfield)

function Base.getindex(f::TransientMultiFieldCellField, index::Integer)
  sub_cellfield = f.cellfield[index]

  sub_derivatives = ()
  for derivative in f.derivatives
    sub_derivative = derivative[index]
    sub_derivatives = (sub_derivatives..., sub_derivative)
  end

  TransientSingleFieldCellField(sub_cellfield, sub_derivatives)
end

function Base.getindex(f::TransientMultiFieldCellField,
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

function Base.iterate(f::TransientMultiFieldCellField)
  iterate(f.transient_single_fields)
end

function Base.iterate(f::TransientMultiFieldCellField, state)
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

####################
# TransientFEBasis #
####################
struct TransientFEBasis{A} <: FEBasis
  febasis::A
  derivatives::Tuple{Vararg{A}}
end

# FEBasis interface
CellData.get_data(f::TransientFEBasis) = get_data(f.febasis)

CellData.get_triangulation(f::TransientFEBasis) = get_triangulation(f.febasis)

CellData.DomainStyle(::Type{<:TransientFEBasis{A}}) where {A} = DomainStyle(A)

function CellData.change_domain(f::TransientFEBasis, trian::Triangulation,
  target_domain::DomainStyle)
  change_domain(f.febasis, trian, target_domain)
end

Fields.gradient(f::TransientFEBasis) = gradient(f.febasis)

Fields.∇∇(f::TransientFEBasis) = ∇∇(f.febasis)

FESpaces.BasisStyle(::Type{<:TransientFEBasis{A}}) where {A} = BasisStyle(A)

# Time derivatives
function ∂t(f::Union{TransientCellField,TransientFEBasis})
  cellfield, derivatives = first_and_tail(f.derivatives)
  TransientCellField(cellfield, derivatives)
end

function ∂tt(f::Union{TransientCellField,TransientFEBasis})
  ∂t(∂t(f))
end

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
