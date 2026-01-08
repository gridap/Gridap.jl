
"""
    from_json(::Type{T},s::AbstractString;kwargs...) where T

Accepts keyword arguments passed to `JSON.parse`.
"""
function from_json(::Type{T},s::AbstractString;kwargs...) where T
  json_dict = JSON.parse(s;kwargs...)
  dict = _decode_json_dict(json_dict)
  from_dict(T,dict)
end

"""
    from_json_file(::Type{T},s::AbstractString;kwargs...) where T

Accepts keyword arguments passed to `JSON.parsefile`.
"""
function from_json_file(::Type{T},s::AbstractString;kwargs...) where T
  json_dict = JSON.parsefile(s;kwargs...)
  dict = _decode_json_dict(json_dict)
  from_dict(T,dict)
end

"""
    to_json(object;kwargs...)

Accepts keyword arguments passed to `JSON.json`.
"""
function to_json(object;kwargs...)
  dict = to_dict(object)
  JSON.json(dict;kwargs...)
end

"""
    to_json_file(object,filename;kwargs...)

Accepts keyword arguments passed to `JSON.json`.
"""
function to_json_file(object,filename;kwargs...)
  dict = to_dict(object)
  open(filename,"w") do f
    JSON.json(f,dict;kwargs...)
  end
end

JSON.lower(object::GridapType) = to_dict(object)

const JSONObject = Union{Dict{String,Any},JSON.Object}

function _decode_json_dict(json_dict::JSONObject)
  _decode_item(x) = x
  _decode_item(x::JSONObject) = _decode_json_dict(x)

  dict = Dict{Symbol,Any}()
  for k in keys(json_dict)
    v = json_dict[k]
    if isa(v,AbstractVector)
      w = collect((_decode_item(vi) for vi in v))
    else
      w = _decode_item(v)
    end
    dict[Symbol(k)] = w
  end
  dict
end
