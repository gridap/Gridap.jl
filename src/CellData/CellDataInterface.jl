
"""
Trait that signals if a CellData type is implemented in the physical or the reference domain
"""
abstract type DomainStyle end
struct ReferenceDomain <: DomainStyle end
struct PhysicalDomain <: DomainStyle end
DomainStyle(::T) where T = DomainStyle(T)

"""
Data associated with the cells of a Triangulation.
CellData objects behave as if they are defined in the physical space of the triangulation.
But in some cases they are implemented as reference quantities plus some transformation to the physical domain.
"""
abstract type CellData <: GridapType end

"""
Get the stored array of cell-wise data. It can be defined in the physical or the reference domain.
"""
get_cell_data(a::CellData) = @abstractmethod

"""
Tell if the stored array is in the reference or physical domain
"""
DomainStyle(::Type{<:CellData}) = @abstractmethod

"""
Return the underlying Triangulation object
"""
get_triangulation(a::CellData) = @abstractmethod

"""
Change the underlying data to the target domain
"""
change_domain(a::CellField,target_domain::DomainStyle) = change_domain(a,DomainStyle(a),target_domain)
change_domain(a::CellField,input_domain::T,target_domain::T) where T<: DomainStyle = a
change_domain(a::CellField,input_domain::DomainStyle,target_domain::DomainStyle) = @abstractmethod

# Tester
"""
"""
function test_cell_data(a::CellData)
  @test isa(get_cell_data(a),AbstractArray)
  @test isa(get_triangulation(a),Triangulation)
  @test isa(DomainStyle(a),DomainStyle)
end

# Some API

"""
Get the raw array of cell data defined in the physical space.
"""
get_array(a::CellData) = get_cell_data(change_domain(a,PhysicalDomain()))

"""
"""
Base.length(a::CellData) = length(get_cell_data(a))


