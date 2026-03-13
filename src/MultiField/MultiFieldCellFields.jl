
struct MultiFieldCellField{DS<:DomainStyle} <: CellField
  single_fields::Vector{<:CellField}
  domain_style::DS

  function MultiFieldCellField(single_fields::Vector{<:CellField})
    @assert !isempty(single_fields)
    is_ref = any(sf -> isa(DomainStyle(sf), ReferenceDomain), single_fields)
    domain_style = ifelse(is_ref, ReferenceDomain(), PhysicalDomain())
    new{typeof(domain_style)}(single_fields,domain_style)
  end
end

function CellData.get_data(f::MultiFieldCellField)
  s = """
  Function get_data is not implemented for MultiFieldCellField at this moment.
  You need to extract the individual fields and then evaluate them separately.

  If ever implemented, evaluating a `MultiFieldCellField` directly would provide,
  at each evaluation point, a tuple with the value of the different fields.
  """
  @notimplemented s
end

function CellData.get_triangulation(f::MultiFieldCellField)
  trian = get_triangulation(first(f.single_fields))
  @check all(sf -> trian === get_triangulation(sf), f.single_fields)
  trian
end

CellData.DomainStyle(::Type{MultiFieldCellField{DS}}) where DS = DS()
num_fields(a::MultiFieldCellField) = length(a.single_fields)
Base.getindex(a::MultiFieldCellField,i::Integer) = a.single_fields[i]
Base.iterate(a::MultiFieldCellField) = iterate(a.single_fields)
Base.iterate(a::MultiFieldCellField,state) = iterate(a.single_fields,state)
Base.length(a::MultiFieldCellField) = num_fields(a)

function LinearAlgebra.dot(a::MultiFieldCellField,b::MultiFieldCellField)
  @check num_fields(a) == num_fields(b)
  return sum(map(dot,a.single_fields,b.single_fields))
end

function Base.show(io::IO,::MIME"text/plain",f::MultiFieldCellField)
  show(io,f)
  print(io,":")
  show_multifield(io,f)
end

function show_multifield(io,f)
  trian = get_triangulation(first(f))
  print(io,"\n num_fields: $(num_fields(f))")
  if all(fi -> trian === get_triangulation(fi), f)
    print(io,"\n num_cells: $(num_cells(trian))")
    print(io,"\n DomainStyle: $(DomainStyle(f))")
    print(io,"\n Triangulation: $(trian)")
    print(io,"\n Triangulation id: $(objectid(trian))")
  else
    print(io,"\n DomainStyle: $(DomainStyle(f))")
    print(io,"\n Triangulations: ")
    for fi in f
      ti = get_triangulation(fi)
      print(io,"\n  > $(ti): num_cells = $(num_cells(ti)), id = $(objectid(ti))")
    end
  end
end

struct MultiFieldFEBasisComponent{B} <: FEBasis
  cell_basis::AbstractArray
  single_field::B
  fieldid::Int
  nfields::Int
end

function block_fields(cell_fields,::TestBasis,fieldid::Integer,nfields::Integer)
  lazy_map(BlockMap(nfields,fieldid),cell_fields)
end
function block_fields(cell_fields,::TrialBasis,fieldid::Integer,nfields::Integer)
  lazy_map(BlockMap((1,nfields),fieldid),cell_fields)
end

function MultiFieldFEBasisComponent(
  single_field::SingleFieldFEBasis,fieldid::Integer,nfields::Integer
)
  cell_basis = block_fields(get_data(single_field),BasisStyle(single_field),fieldid,nfields)
  MultiFieldFEBasisComponent(cell_basis,single_field,fieldid,nfields)
end

function MultiFieldFEBasisComponent(
  single_field::CellField,fieldid::Integer,nfields::Integer
)
  function basis_style(cf)
    FT = eltype(CellData.get_data(cf))
    (FT <: AbstractVector{<:Field}) && (return TestBasis())
    (FT <: AbstractMatrix{<:Field}) && (return TrialBasis())
    error("Cannot determine BasisStyle from CellField.")
  end
  cell_basis = block_fields(get_data(single_field),basis_style(single_field),fieldid,nfields)
  MultiFieldFEBasisComponent(cell_basis,single_field,fieldid,nfields)
end

function MultiFieldFEBasisComponent(
  single_field::CellFieldAt{S,<:SingleFieldFEBasis}, fieldid::Integer, nfields::Integer
) where S
  sf = single_field.parent
  mf = MultiFieldFEBasisComponent(sf,fieldid,nfields)
  CellFieldAt{S}(mf)
end

function MultiFieldFEBasis(
  single_fields::Vector{<:CellField},
  field_ids = 1:length(single_fields),
  nfields = length(single_fields),
)
  @assert length(single_fields) == length(field_ids)
  blocked_fields = map(single_fields, field_ids) do sf, fid
    MultiFieldFEBasisComponent(sf,fid,nfields)
  end
  MultiFieldCellField(blocked_fields)
end

# Swapping field ids

function MultiFieldFEBasis(cf::MultiFieldCellField, args...)
  MultiFieldFEBasis(cf.single_fields, args...)
end

function MultiFieldFEBasis(cf::Tuple, args...)
  MultiFieldFEBasis([cf...], args...)
end

function MultiFieldFEBasisComponent(
  single_field::MultiFieldFEBasisComponent,fieldid::Integer,nfields::Integer
)
  MultiFieldFEBasisComponent(single_field.single_field,fieldid,nfields)
end

# CellField API 

CellData.get_data(f::MultiFieldFEBasisComponent) = f.cell_basis
CellData.get_triangulation(f::MultiFieldFEBasisComponent) = get_triangulation(f.single_field)
FESpaces.BasisStyle(::Type{<:MultiFieldFEBasisComponent{B}}) where B = BasisStyle(B)
CellData.DomainStyle(::Type{<:MultiFieldFEBasisComponent{B}}) where B = DomainStyle(B)

function FESpaces.CellData.similar_cell_field(
  f::MultiFieldFEBasisComponent,cell_data,trian,ds::DomainStyle)
  @notimplemented
end
function FESpaces.similar_fe_basis(
  f::MultiFieldFEBasisComponent,cell_data,trian,bs::BasisStyle,ds::DomainStyle)
  @notimplemented
end

for fun in (:gradient,:DIV,:∇∇)
  @eval begin
    function $fun(f::MultiFieldFEBasisComponent)
      g = $fun(f.single_field)
      MultiFieldFEBasisComponent(g,f.fieldid,f.nfields)
    end
  end
end

function CellData.change_domain(
  a::MultiFieldFEBasisComponent,tdomain::DomainStyle
)
  sf = change_domain(a.single_field,tdomain)
  MultiFieldFEBasisComponent(sf,a.fieldid,a.nfields)
end

function CellData.change_domain(
  a::MultiFieldFEBasisComponent, ttrian::Triangulation, tdomain::DomainStyle
)
  sf = change_domain(a.single_field,ttrian,tdomain)
  MultiFieldFEBasisComponent(sf,a.fieldid,a.nfields)
end

function CellData.change_domain(
  a::CellFieldAt{S,<:MultiFieldFEBasisComponent},
  tdomain::DomainStyle
) where S
  mf = a.parent
  sfin = CellFieldAt{S}(mf.single_field)
  sfout = change_domain(sfin,tdomain)
  MultiFieldFEBasisComponent(sfout,mf.fieldid,mf.nfields)
end

function CellData.change_domain(
  a::CellFieldAt{S,<:MultiFieldFEBasisComponent},
  ttrian::Triangulation,
  tdomain::DomainStyle
) where S
  mf = a.parent
  sfin = CellFieldAt{S}(mf.single_field)
  sfout = change_domain(sfin,ttrian,tdomain)
  MultiFieldFEBasisComponent(sfout,mf.fieldid,mf.nfields)
end
