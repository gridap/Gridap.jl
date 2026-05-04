
"""
    to_dict(object) -> Dict

Serialize `object` into a dictionary of type `Dict{Symbol,Any}`.

Values stored in this dictionary must be native Julia data types
(e.g., `Real`, `Integer`, `String`, `Vector`, `Dict`).
Dictionary keys are `Symbols`. This is typically the first step
in serializing an object to formats like JSON or BSON.
"""
function to_dict(object)
  @abstractmethod
end

"""
    check_dict(::Type{T}, dict::Dict) where T

Check if the dictionary `dict` is valid for an object of type `T`.

Returns successfully if the dictionary is compatible with the type `T`,
otherwise throws an error.
"""
function check_dict(::Type{T},dict::Dict) where T
  @abstractmethod
end

"""
    from_dict(::Type{T}, dict::Dict) where T

De-serialize an object of type `T` from the dictionary `dict`.

Values in the dictionary should be native Julia data types.
Dictionary keys are `Symbols`.
"""
function from_dict(::Type{T},dict::Dict) where T
  @abstractmethod
end
