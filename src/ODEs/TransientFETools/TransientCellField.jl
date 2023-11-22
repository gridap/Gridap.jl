######################
# TransientCellField #
######################
abstract type TransientCellField <: CellField end

# CellField interface
get_data(f::TransientCellField) = @abstractmethod

get_triangulation(f::TransientCellField) = @abstractmethod

DomainStyle(::Type{TransientCellField}) = @abstractmethod

gradient(f::TransientCellField) = @abstractmethod

∇∇(f::TransientCellField) = @abstractmethod

function change_domain(f::TransientCellField, trian::Triangulation,
  target_domain::DomainStyle)
  @abstractmethod
end

#################################
# TransientSingleFieldCellField #
#################################
struct TransientSingleFieldCellField{A} <: TransientCellField
  cellfield::A
  derivatives::Tuple # {Vararg{A,B} where B}
end

const SingleFieldTypes = Union{GenericCellField,SingleFieldFEFunction}

function TransientCellField(single_field::SingleFieldTypes, derivatives::Tuple)
  TransientSingleFieldCellField(single_field, derivatives)
end

# CellField interface
get_data(f::TransientSingleFieldCellField) = get_data(f.cellfield)

get_triangulation(f::TransientSingleFieldCellField) = get_triangulation(f.cellfield)

DomainStyle(::Type{<:TransientSingleFieldCellField{A}}) where {A} = DomainStyle(A)

gradient(f::TransientSingleFieldCellField) = gradient(f.cellfield)

∇∇(f::TransientSingleFieldCellField) = ∇∇(f.cellfield)

function change_domain(f::TransientSingleFieldCellField, trian::Triangulation,
  target_domain::DomainStyle)
  change_domain(f.cellfield, trian, target_domain)
end

# Skeleton-related operations
function getproperty(f::TransientSingleFieldCellField, sym::Symbol)
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

####################
# TransientFEBasis #
####################
struct TransientFEBasis{A} <: FEBasis
  febasis::A
  derivatives::Tuple{Vararg{A}}
end

# FEBasis methods
get_data(f::TransientFEBasis) = get_data(f.febasis)

get_triangulation(f::TransientFEBasis) = get_triangulation(f.febasis)

DomainStyle(::Type{<:TransientFEBasis{A}}) where {A} = DomainStyle(A)

BasisStyle(::Type{<:TransientFEBasis{A}}) where {A} = BasisStyle(A)

gradient(f::TransientFEBasis) = gradient(f.febasis)

∇∇(f::TransientFEBasis) = ∇∇(f.febasis)

function change_domain(f::TransientFEBasis, trian::Triangulation,
  target_domain::DomainStyle)
  change_domain(f.febasis, trian, target_domain)
end

# Time derivatives
function ∂t(f::Union{TransientCellField,TransientFEBasis})
  cellfield, derivatives = first_and_tail(f.derivatives)
  TransientCellField(cellfield, derivatives)
end

function ∂tt(f::Union{TransientCellField,TransientFEBasis})
  ∂t(∂t(f))
end
