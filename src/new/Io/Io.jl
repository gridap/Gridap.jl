"""
The exported names in this module are:

$(EXPORTS)
"""
module Io

using DocStringExtensions

using Gridap.Helpers
import JSON

export to_dict
export from_dict

export to_json
export to_json_file
export from_json
export from_json_file

include("IoInterfaces.jl")

include("Json.jl")

end
