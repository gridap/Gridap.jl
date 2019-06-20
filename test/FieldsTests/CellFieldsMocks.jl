module CellFieldsMocks

using Gridap
using Gridap.CellFieldsOperations: _merge_val_and_grad

using ..CellValuesMocks
using ..FieldsMocks

export IterCellFieldMock
export IndexCellFieldMock
export IterCellBasisMock
export IndexCellBasisMock
export IterCellGeomapMock
export IndexCellGeomapMock

function IterCellFieldMock(D::Int,X::Type,l::Int)
  f = MockField(D,X)
  g = gradient(f)
  cf = TestIterCellValue(f,l)
  cg = TestIterCellValue(g,l)
  _merge_val_and_grad(cf,cg)
end

function IndexCellFieldMock(D::Int,X::Type,l::Int)
  f = MockField(D,X)
  g = gradient(f)
  cf = TestIndexCellValue(f,l)
  cg = TestIndexCellValue(g,l)
  _merge_val_and_grad(cf,cg)
end

function IterCellBasisMock(D::Int,X::Type,l::Int)
  f = MockBasis(D,X)
  g = gradient(f)
  cf = TestIterCellValue(f,l)
  cg = TestIterCellValue(g,l)
  _merge_val_and_grad(cf,cg)
end

function IndexCellBasisMock(D::Int,X::Type,l::Int)
  f = MockBasis(D,X)
  g = gradient(f)
  cf = TestIndexCellValue(f,l)
  cg = TestIndexCellValue(g,l)
  _merge_val_and_grad(cf,cg)
end

function IterCellGeomapMock(D::Int,X::Type,l::Int)
  f = MockGeomap(D,X)
  g = gradient(f)
  cf = TestIterCellValue(f,l)
  cg = TestIterCellValue(g,l)
  _merge_val_and_grad(cf,cg)
end

function IndexCellGeomapMock(D::Int,X::Type,l::Int)
  f = MockGeomap(D,X)
  g = gradient(f)
  cf = TestIndexCellValue(f,l)
  cg = TestIndexCellValue(g,l)
  _merge_val_and_grad(cf,cg)
end

end # module
