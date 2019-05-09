module CellMaps

using Reexport

include("AbstractCellMaps.jl")

@reexport using Numa.CellMaps.AbstractCellMaps

include("CellMapValues.jl")

include("ConstantCellMaps.jl")

include("Testers.jl")

end #module CellMaps

#using Numa.Helpers
#
#using Numa.Maps
#using Numa.Maps: MapFromUnaryOp
#using Numa.Maps: MapFromBinaryOp
#using Numa.Maps: FieldFromExpand
#using Numa.Maps: FieldFromCompose
#using Numa.Maps: FieldFromComposeExtended
#
#using Numa.FieldValues
#import Numa.FieldValues: inner, outer
#using Numa.CellValues
#using Numa.CellValues: CellArrayFromUnaryOp
#using Numa.CellValues: CellArrayFromBroadcastUnaryOp
#using Numa.CellValues: CellArrayFromBroadcastBinaryOp
#using Numa.CellValues: CachedArray
#using Numa.CellValues: setsize!
#
#export CellMap
#export CellField
#export CellBasis
#export CellGeomap
#
#export CellMapValues
#export CellBasisValues
#export CellPoints
#
#export ConstantCellMap
#
#export expand
#export varinner
#export attachgeomap
#export compose
#
#import Base: +, -, *, /, ∘
#import Base: iterate
#import Base: length
#import Base: eltype
#import Base: size
#import Base: getindex, setindex!
#
#import Numa: evaluate, gradient, ∇
#import Numa: evaluate!, return_size
#import Numa.CellValues: cellsize
#import Numa.CellValues: inputcellarray, computesize, computevals!

#include("ConcreteCellMaps.jl")
#include("Operators.jl")
#include("CellBasisWithGeomap.jl")
