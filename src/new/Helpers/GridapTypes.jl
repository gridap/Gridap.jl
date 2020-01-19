
"""
    abstract type GridapType end
"""
abstract type GridapType end

function Base.show(io::IO,object::GridapType)
  print(io,"$(nameof(typeof(object)))()")
end

