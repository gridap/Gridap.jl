
"""
    abstract type DomainStyle

Trait that signals if a CellDatum type is implemented in the physical or the
reference domain, the possible values are `ReferenceDomain()` and
`PhysicalDomain()`.
"""
abstract type DomainStyle end
struct ReferenceDomain <: DomainStyle end
struct PhysicalDomain <: DomainStyle end

"""
    abstract type CellDatum <: GridapType

Data associated with the cells of a Triangulation.
CellDatum objects behave as if they are defined in the physical space of the triangulation.
But in some cases they are implemented as reference quantities plus some transformation to the physical domain.
"""
abstract type CellDatum <: GridapType end

"""
    get_data(a::CellDatum)

Get the stored array of cell-wise data. It can be defined in the physical or the reference domain.
"""
get_data(a::CellDatum) = @abstractmethod

"""
    DomainStyle(::Type{<:CellDatum})

Tell if the stored array is in the reference or physical domain.
"""
DomainStyle(::Type{<:CellDatum}) = @abstractmethod

"""
    get_triangulation(a::CellDatum)

Return the underlying Triangulation object.
"""
get_triangulation(a::CellDatum) = @abstractmethod

"""
    change_domain(a::CellDatum, target_domain)
    change_domain(a::CellDatum, input_domain, target_domain)

where `a` isa [`CellDatum`](@ref) and the domains are [`DomainStyle`](@ref)s.
Change the underlying data to the target domain
"""
change_domain(a::CellDatum,target_domain::DomainStyle) = change_domain(a,DomainStyle(a),target_domain)
change_domain(a::CellDatum,input_domain::T,target_domain::T) where T<: DomainStyle = a
change_domain(a::CellDatum,input_domain::DomainStyle,target_domain::DomainStyle) = @abstractmethod

# Tester
"""
"""
function test_cell_datum(a::CellDatum)
  @test isa(get_data(a),AbstractArray)
  @test isa(get_triangulation(a),Triangulation)
  @test isa(DomainStyle(a),DomainStyle)
end

# Some API

DomainStyle(::T) where T<:CellDatum = DomainStyle(T)

"""
    get_array(a::CellDatum)

Get the raw array of `a` datas defined in the physical space.
"""
get_array(a::CellDatum) = get_data(change_domain(a,PhysicalDomain()))

"""
    num_cells(a::CellDatum)

Number of cells of `a`.
"""
num_cells(a::CellDatum) = num_cells(get_triangulation(a))
