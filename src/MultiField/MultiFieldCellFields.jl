
struct MultiFieldCellField{DS<:DomainStyle} <: CellField
  single_fields::Vector{<:CellField}
  trian::Triangulation
  domain_style::DS

  function MultiFieldCellField(single_fields::Vector{<:CellField})
    @assert length(single_fields) > 0
    f1 = first(single_fields)

    # Find a suitable domain style
    if any( map(i->DomainStyle(i)==ReferenceDomain(),single_fields) )
      domain_style = ReferenceDomain()
    else
      domain_style = PhysicalDomain()
    end

    trian = get_triangulation(f1)
    #@check all(map(i-> get_triangulation(i) === trian,single_fields))

    new{typeof(domain_style)}(single_fields,trian,domain_style)
  end
end

function CellData.get_data(f::MultiFieldCellField)
  s = """
  Function get_data is not implemented for MultiFieldCellField at this moment.
  You need to extract the individual fields and then evaluate them separatelly.

  If ever implement this, evaluating a `MultiFieldCellField` directly would provide,
  at each evaluation point, a tuple with the value of the different fields.
  """
  @notimplemented s
end

CellData.get_triangulation(f::MultiFieldCellField) = f.trian
CellData.DomainStyle(::Type{MultiFieldCellField{DS}}) where DS = DS()
num_fields(a::MultiFieldCellField) = length(a.single_fields)
Base.getindex(a::MultiFieldCellField,i::Integer) = a.single_fields[i]
Base.iterate(a::MultiFieldCellField)  = iterate(a.single_fields)
Base.iterate(a::MultiFieldCellField,state)  = iterate(a.single_fields,state)
