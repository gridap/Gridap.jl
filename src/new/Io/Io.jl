"""
The exported names in this module are:

$(EXPORTS)
"""
module Io

using DocStringExtensions
using JLD2
using FileIO

using Gridap.Helpers
import JSON
using BSON

export to_dict
export from_dict

export to_json
export to_json_file
export to_bson_file
export to_jld2_file
export from_json
export from_json_file
export from_bson_file
export from_jld2_file

include("IoInterfaces.jl")

include("Json.jl")

include("Bson.jl")

include("JLD2.jl")
end
