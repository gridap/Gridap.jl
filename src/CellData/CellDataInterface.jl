
"""
Trait that signals if a CellDatum type is implemented in the physical or the reference domain
"""
abstract type DomainStyle end
struct ReferenceDomain <: DomainStyle end
struct PhysicalDomain <: DomainStyle end

"""
Data associated with the cells of a Triangulation.
CellDatum objects behave as if they are defined in the physical space of the triangulation.
But in some cases they are implemented as reference quantities plus some transformation to the physical domain.
"""
abstract type CellDatum <: GridapType end

"""
Get the stored array of cell-wise data. It can be defined in the physical or the reference domain.
"""
get_data(a::CellDatum) = @abstractmethod

"""
Tell if the stored array is in the reference or physical domain
"""
DomainStyle(::Type{<:CellDatum}) = @abstractmethod

"""
Return the underlying Triangulation object
"""
get_triangulation(a::CellDatum) = @abstractmethod

"""
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
Get the raw array of cell data defined in the physical space.
"""
get_array(a::CellDatum) = get_data(change_domain(a,PhysicalDomain()))

"""
"""
#Base.length(a::CellDatum) = num_cells(get_triangulation(a))
num_cells(a::CellDatum) = num_cells(get_triangulation(a))
