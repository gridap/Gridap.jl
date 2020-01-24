module FieldsBenchs

include("MockFieldsBenchs.jl")

include("ConstantFieldsBenchs.jl")

include("FieldApplyBenchs.jl")

include("FieldArraysBenchs.jl")

include("LincombBenchs.jl")

include("ComposeBenchs.jl")

include("FieldOperationsBenchs.jl")

include("AttachmapBenchs.jl")

include("IntegrateBenchs.jl")

end # module
