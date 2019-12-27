module Io

using Gridap.Helpers
import JSON

export to_dict
export from_dict

export to_json
export to_json_file
export from_json
export from_json_file
export decode_json_dict

include("IoInterfaces.jl")

include("Json.jl")

end
