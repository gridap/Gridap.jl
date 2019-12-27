
"""
    decode_json_dict(::Type{T},json_dict::Dict{String,Any}) where T

Tranform a dictionary resulting from the parsing of a JSON string/file into another
one such that

    dict = decode_json_dict(T,json_dict)
    object = from_dict(T, dict)

Returns an object of type `T`.
"""
function decode_json_dict(::Type{T},dict::Dict{String,Any}) where T
  @abstractmethod
end

function from_json(::Type{T},s::AbstractString) where T
  json_dict = JSON.parse(s)
  dict = decode_json_dict(T,json_dict)
  from_dict(T,dict)
end

function from_json_file(::Type{T},s::AbstractString) where T
  json_dict = JSON.parsefile(s)
  dict = decode_json_dict(T,json_dict)
  from_dict(T,dict)
end

"""
    to_json(object)
"""
function to_json(object)
  dict = to_dict(object)
  JSON.json(dict)
end

"""
    to_json_file(object,filename)
"""
function to_json_file(object,filename)
  dict = to_dict(object)
  open(filename,"w") do f
    JSON.print(f,dict)
  end
end

JSON.lower(object::GridapType) = to_dict(object)
