"""
    from_bson_file(::Type{T}, path::AbstractString) where T

De-serializes a BSON file into an object of type `T`.
"""
function from_bson_file(::Type{T},s::AbstractString) where T
  dict = BSON.load(s)
  from_dict(T,dict)
end

"""
    to_bson_file(object, path::AbstractString)

Serializes `object` into a BSON file.
"""
function to_bson_file(object,filename)
  dict = to_dict(object)
  bson(filename, dict)
end

