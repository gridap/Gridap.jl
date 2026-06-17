module ReferenceFEs

using Gridap.Helpers
using Gridap.TensorValues
using Gridap.Polynomials
using Gridap.ReferenceFEs

using Gridap.ReferenceFEs: ReferenceFEName
using Gridap.ReferenceFEs: Conformity
using Gridap.ReferenceFEs: GenericRefFE

import Gridap.ReferenceFEs: ReferenceFE
import Gridap.ReferenceFEs: valid_conformity_symbols

export Loop
export LoopRefFE
export LoopConformity
export loop

include("LoopRefFEs.jl")

end
