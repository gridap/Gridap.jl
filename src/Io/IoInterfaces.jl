
"""
    to_dict(object) -> Dict

Serialize `object` into a dictionary of type `Dict{Symbol,Any}`.
Values stored into this dictionary must be of any native 
Julia data type (Real, Integer, String, etc.)
Dictionary keys are `Symbols`.
"""
function to_dict(object)
  @abstractmethod
end

"""
    check_dict(::Type{T},dict::Dict) where T

Check validity of a dictionary `dict` for an object of type `T`.
It runs successfully if the dictionary is valid for a particular 
type or throws an error in any other case.
"""
function check_dict(::Type{T},dict::Dict) where T
  @abstractmethod
end

"""
    from_dict(::Type{T},dict::Dict) where T

De-serialize an object of type `T` from the dictionary `dict`.
Values stored into this dictionary must be of any native 
Julia data type (Real, Integer, String, etc.)
Dictionary keys are `Symbols`.
"""
function from_dict(::Type{T},dict::Dict) where T
  @abstractmethod
end
