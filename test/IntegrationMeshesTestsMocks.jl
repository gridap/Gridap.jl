
# TODO: the following is only a dummy implementation in order to
# test the interfaces at this early stage.
# It will be removed in the future and replaced by other more general
# and elegant implementations

include("PolynomialsTestsMocks.jl")

using Numa.CellArrays: IndexableCellArray
import Numa.CellArrays: cellsize
using Numa.IntegrationMeshes
import Numa.IntegrationMeshes: cellcoordinates, cellbasis

struct DummyCellCoordinates2D <: IndexableCellArray{Point{2},1}
  x::Array{Point{2},2}
  c::Array{Point{2},1}
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
  c = Array{Point{2},1}(undef,(4,))
  DummyCellCoordinates2D(x,c)
end

function Base.getindex(self::DummyCellCoordinates2D,cell::Int)
  @inbounds for i in 1:4
    self.c[i] = self.x[i,cell]
  end
  self.c
end

Base.length(self::DummyCellCoordinates2D) = size(self.x,2)

cellsize(self::DummyCellCoordinates2D) = (4,)

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
