######################
# TransientCellField #
######################
"""
    abstract type TransientCellField <: CellField end

Transient version of `CellField`.

# Mandatory
- [`time_derivative(f)`](@ref)
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

# TransientCellField interface
function time_derivative(f::TransientCellField)
  @abstractmethod
end

#################################
# TransientSingleFieldCellField #
#################################
"""
    struct TransientSingleFieldCellField <: TransientCellField end

Transient `CellField` for a single-field `FESpace`.
"""
struct TransientSingleFieldCellField{A} <: TransientCellField
  cellfield::A
  derivatives::Tuple # {Vararg{A,B} where B}
end

# Default constructor (see `TransientMultiFieldCellField` for the implementations
# of `TransientCellField` when the field is a `MultiFieldCellField`)
function TransientCellField(field::CellField, derivatives::Tuple)
  TransientSingleFieldCellField(field, derivatives)
end

# CellField interface
CellData.get_data(f::TransientSingleFieldCellField) = get_data(f.cellfield)

CellData.get_triangulation(f::TransientSingleFieldCellField) = get_triangulation(f.cellfield)

CellData.DomainStyle(::Type{<:TransientSingleFieldCellField{A}}) where {A} = DomainStyle(A)

function CellData.change_domain(
  f::TransientSingleFieldCellField, trian::Triangulation,
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

# TransientCellField interface
function time_derivative(f::TransientSingleFieldCellField)
  cellfield, derivatives = first_and_tail(f.derivatives)
  TransientCellField(cellfield, derivatives)
end

################################
# TransientMultiFieldCellField #
################################
"""
    struct TransientMultiFieldCellField <: TransientCellField end

Transient `CellField` for a multi-field `FESpace`.
"""
struct TransientMultiFieldCellField{A} <: TransientCellField
  cellfield::A
  derivatives::Tuple
  transient_single_fields::Vector{<:TransientCellField} # used to iterate
end

const MultiFieldTypes = Union{MultiFieldCellField,MultiFieldFEFunction}
function TransientMultiFieldCellField(fields::MultiFieldTypes, derivatives::Tuple)
  _flat = _to_transient_single_fields(fields, derivatives)
  TransientMultiFieldCellField(fields, derivatives, _flat)
end

# Default constructors
function TransientCellField(fields::MultiFieldTypes, derivatives::Tuple)
  TransientMultiFieldCellField(fields, derivatives)
end

function TransientCellField(fields::TransientMultiFieldCellField, derivatives::Tuple)
  TransientMultiFieldCellField(fields, derivatives)
end

# CellField interface
function CellData.get_data(f::TransientMultiFieldCellField)
  s = """
  Function `get_data` is not implemented for `TransientMultiFieldCellField` at
  this moment. You need to extract the individual fields and then evaluate them
  separately.

  If this function is ever to be implemented, evaluating a `MultiFieldCellField`
  directly would provide, at each evaluation point, a tuple with the value of
  the different fields.
  """
  @notimplemented s
end

CellData.get_triangulation(f::TransientMultiFieldCellField) = get_triangulation(f.cellfield)

CellData.DomainStyle(::Type{TransientMultiFieldCellField{A}}) where {A} = DomainStyle(A)

function CellData.change_domain(
  f::TransientMultiFieldCellField, trian::Triangulation,
  target_domain::DomainStyle
)
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

function Base.getindex(
  f::TransientMultiFieldCellField,
  indices::AbstractVector{<:Integer}
)
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

# TransientCellField interface
function time_derivative(f::TransientMultiFieldCellField)
  cellfield, derivatives = first_and_tail(f.derivatives)

  single_field_derivatives = map(cellfield, derivatives...) do cellfield, derivatives...
    TransientSingleFieldCellField(cellfield, derivatives)
  end

  TransientMultiFieldCellField(
    cellfield, derivatives,
    single_field_derivatives
  )
end

####################
# TransientFEBasis #
####################
"""
    struct TransientFEBasis <: FEBasis end

Transient `FEBasis`.
"""
struct TransientFEBasis{A} <: FEBasis
  febasis::A
  derivatives::Tuple{Vararg{A}}
end

# CellField interface
CellData.get_data(f::TransientFEBasis) = get_data(f.febasis)

CellData.get_triangulation(f::TransientFEBasis) = get_triangulation(f.febasis)

CellData.DomainStyle(::Type{<:TransientFEBasis{A}}) where {A} = DomainStyle(A)

function CellData.change_domain(
  f::TransientFEBasis, trian::Triangulation,
  target_domain::DomainStyle
)
  change_domain(f.febasis, trian, target_domain)
end

Fields.gradient(f::TransientFEBasis) = gradient(f.febasis)

Fields.∇∇(f::TransientFEBasis) = ∇∇(f.febasis)

# FEBasis interface
FESpaces.BasisStyle(::Type{<:TransientFEBasis{A}}) where {A} = BasisStyle(A)

# Transient FEBasis interface
function time_derivative(f::TransientFEBasis)
  cellfield, derivatives = first_and_tail(f.derivatives)
  TransientCellField(cellfield, derivatives)
end

#########
# Utils #
#########
"""
    _to_transient_single_fields(
      multi_field,
      derivatives
    ) -> Vector{<:TransientSingleFieldCellField}

Convert a `TransientMultiFieldCellField` into a vector of
`TransientSingleFieldCellField`s.
"""
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
