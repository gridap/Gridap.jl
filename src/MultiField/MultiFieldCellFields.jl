struct MultiFieldCellField{DS<:DomainStyle} <: CellField
  single_fields::Vector{<:CellField}
  domain_style::DS

  function MultiFieldCellField(single_fields::Vector{<:CellField})
    @assert length(single_fields) > 0
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
