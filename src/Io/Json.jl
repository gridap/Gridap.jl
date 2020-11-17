
"""
    from_json(::Type{T},s::AbstractString) where T
"""
function from_json(::Type{T},s::AbstractString) where T
  json_dict = JSON.parse(s;inttype=Int)
  dict = _decode_json_dict(json_dict)
  from_dict(T,dict)
end

"""
    from_json_file(::Type{T},s::AbstractString) where T
"""
function from_json_file(::Type{T},s::AbstractString) where T
  json_dict = JSON.parsefile(s;inttype=Int)
  dict = _decode_json_dict(json_dict)
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

function _decode_json_dict(json_dict::Dict{String,Any})
  dict = Dict{Symbol,Any}()
  for k in keys(json_dict)
    v = json_dict[k]
    if isa(v,Dict{String,Any})
      w = _decode_json_dict(v)
    elseif isa(v,AbstractVector)
      _w = []
      for vi in v
        if isa(vi,Dict{String,Any})
          wi = _decode_json_dict(vi)
        else
          wi = vi
        end
        push!(_w,wi)
      end
      w = collect(_w)
    else
      w = v
    end
    dict[Symbol(k)] = w
  end
  dict
end

