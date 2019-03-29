
"""
An array of values for each cell
"""
const CellValues{T} = CellArray{T,1} where T

"""
Abstract type that represents a field with value of type T
evaluated at a collection of points in each cell
"""
const CellFieldValues{T} = CellValues{T} where T <: FieldValue

"""
Abstract type that represents a function basis with value of type T
evaluated at a collection of points in each cell
"""
const CellBasisValues{T} = CellArray{T,2} where T <: FieldValue

"""
An array of points for each cell
"""
const CellPoints{D} = CellFieldValues{Point{D}} where D

# Operations

inner(a::CellFieldValues,b::CellFieldValues) = binner(a,b)

inner(a::CellBasisValues,b::CellFieldValues) = binner(a,cellnewaxis(b,dim=1))

inner(a::CellBasisValues,b::CellBasisValues) = binner(cellnewaxis(a,dim=2),cellnewaxis(b,dim=1))

expand(a::CellBasisValues,b::CellFieldValues) = cellsum(bouter(a,cellnewaxis(b,dim=2)),dim=1)

