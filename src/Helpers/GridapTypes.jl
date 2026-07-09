
"""
    abstract type GridapType end
"""
abstract type GridapType end

function Base.show(io::IO,object::GridapType)
  print(io,"$(nameof(typeof(object)))()")
end

# Integer types

"""
    const GridapLocalInt <: Integer

Integer type used for local indexing in Gridap. It will be used, for instance, for the 
local numbering of faces and degrees of freedom within a cell. 

Defaults to `Int8`, which is enough for most applications (up to 127 vertices/faces/dofs per cell).
Its value can be changed at compile time using [`set_local_integer_type`](@ref).
"""
const GridapLocalInt = eval(Symbol(@load_preference("local_integer_type", "Int8")))

"""
    set_local_integer_type(T::Type{<:Integer})

Sets the integer type used for local indexing in Gridap. See [`GridapLocalInt`](@ref) for more details.
"""
function set_local_integer_type(T::Type{<:Integer})
  @set_preferences!("local_integer_type" => string(T))
  @info("New local integer type set; restart your Julia session for this change to take effect!")
end
