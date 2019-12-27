
"""
    to_dict(object) -> Dict

Serialize `object` into a dictionary
"""
function to_dict(object)
  @abstractmethod
end

"""
    from_dict(::Type{T},dict::Dict) where T

De-serialize an object of type `T` from the dictionary `dict`.
"""
function from_dict(::Type{T},dict::Dict) where T
  @abstractmethod
end
