export IntegrationDomain, IntegrationMesh, geomap, cellcoordinates, cellbasis, ncells

"""
This is the very minimum needed to describe the
domain for numerical integration
"""
abstract type IntegrationDomain{Z,D} end

geomap(::IntegrationDomain{Z,D} where {Z,D})::CellField{Z,Point{D}}= @abstractmethod

"""
Minimal interface for a mesh used for numerical integration
"""
abstract type IntegrationMesh{Z,D} <: IntegrationDomain{Z,D} end

cellcoordinates(::IntegrationMesh{Z,D} where {Z,D})::CellPoints{D} = @abstractmethod

cellbasis(::IntegrationMesh{Z,D} where {Z,D})::CellBasis{Float64,Z} = @abstractmethod

function geomap(self::IntegrationMesh{Z,D}) where {Z,D}
  coords = cellcoordinates(self)
  basis = cellbasis(self)
  CellFieldFromInterpolation{Z,Point{D},Float64,Point{D}}(basis,coords)
end

function ncells(self::IntegrationMesh)
  coords = cellcoordinates(self)
  length(coords)
end

# Concrete implementations

# TODO: the following is only a dummy implementation in order to
# test the interfaces at this early stage.
# It will be removed in the future and replaced by other more general
# and elegant implementations

struct DummyCellCoordinates2D <: IndexableCellArray{Point{2},1}
  x::Array{Point{2},2}
end

function DummyCellCoordinates2D(;partition::Tuple{Int,Int})
  ncells = prod(partition)
  x = Array{Point{2},2}(undef,(4,ncells))
  pc = MPoint{2}(zeros(2))
  pn = MPoint{2}(zeros(2))
  hx = 1.0/partition[1]
  hy = 1.0/partition[2]
  cell = 1
  for cy in 1:partition[2]
    pc[2] = (cy-1)*hy
    for cx in 1:partition[1]
      pc[1] = (cx-1)*hx
      node = 1
      for ny in 1:2
        pn[2] = pc[2] + (ny-1)*hy
        for nx in 1:2
          pn[1] = pc[1] + (nx-1)*hx
          x[node,cell] = pn
          node += 1
        end
      end
      cell += 1
    end
  end
  DummyCellCoordinates2D(x)
end

Base.getindex(self::DummyCellCoordinates2D,cell::Int) = @view self.x[:,cell]

Base.length(self::DummyCellCoordinates2D) = size(self.x,2)

maxsize(self::DummyCellCoordinates2D) = (4,)

struct DummyIntegrationMesh2D <: IntegrationMesh{2,2}
  cellcoords::DummyCellCoordinates2D
  cellbasis::CellBasisFromSingleInterpolation{2,Float64}
end

function DummyIntegrationMesh2D(;partition::Tuple{Int,Int})
  basis = ShapeFunctionsScalarQua4()
  cellbasis = CellBasisFromSingleInterpolation(basis)
  cellcoords = DummyCellCoordinates2D(partition=partition)
  DummyIntegrationMesh2D(cellcoords,cellbasis)
end

cellcoordinates(self::DummyIntegrationMesh2D) = self.cellcoords

cellbasis(self::DummyIntegrationMesh2D) = self.cellbasis


