module CompressedCellFields

using Gridap
import Gridap: gradient

const CompressedCellFieldLike{D,T,N,R<:FieldLike{D,T,N}} = CompressedCellValue{R}

function gradient(f::CompressedCellFieldLike)
  broadcast(gradient,f)
end

end # module
