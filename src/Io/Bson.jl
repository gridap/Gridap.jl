"""
    function from_bson_file(::Type{T},s::AbstractString) where T
"""
function from_bson_file(::Type{T},s::AbstractString) where T
  dict = BSON.load(s)
  from_dict(T,dict)
end

"""
    to_bson_file(object,filename)
"""
function to_bson_file(object,filename)
  dict = to_dict(object)
  bson(filename, dict)
end

