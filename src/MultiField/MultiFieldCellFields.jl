"""
    struct MultiFieldCellField{DS<:DomainStyle} <: CellField
"""
struct MultiFieldCellField{DS<:DomainStyle} <: CellField
  single_fields::Vector{<:CellField}
  domain_style::DS

  """
      MultiFieldCellField(single_fields::Vector{<:CellField})
  """
  function MultiFieldCellField(single_fields::Vector{<:CellField})
    @assert length(single_fields) > 0
    f1 = first(single_fields)

    # Find a suitable domain style
    if any( map(i->DomainStyle(i)==ReferenceDomain(),single_fields) )
      domain_style = ReferenceDomain()
    else
      domain_style = PhysicalDomain()
    end

    new{typeof(domain_style)}(single_fields,domain_style)
  end
end

function CellData.get_data(f::MultiFieldCellField)
  s = """
  Function get_data is not implemented for MultiFieldCellField at this moment.
  You need to extract the individual fields and then evaluate them separately.

  If ever implement this, evaluating a [`MultiFieldCellField`](@ref) directly
  would provide, at each evaluation point, a tuple with the value of the
  different fields.
  """
  @notimplemented s
end

function CellData.get_triangulation(f::MultiFieldCellField)
  s1 = first(f.single_fields)
  trian = get_triangulation(s1)
  @check all(map(i->trian===get_triangulation(i),f.single_fields))
  trian
end
CellData.DomainStyle(::Type{MultiFieldCellField{DS}}) where DS = DS()
num_fields(a::MultiFieldCellField) = length(a.single_fields)
Base.getindex(a::MultiFieldCellField,i::Integer) = a.single_fields[i]
Base.iterate(a::MultiFieldCellField)  = iterate(a.single_fields)
Base.iterate(a::MultiFieldCellField,state)  = iterate(a.single_fields,state)
Base.length(a::MultiFieldCellField) = num_fields(a)

function LinearAlgebra.dot(a::MultiFieldCellField,b::MultiFieldCellField)
  @check num_fields(a) == num_fields(b)
  return sum(map(dot,a.single_fields,b.single_fields))
end
