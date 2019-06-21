module ConstantCellFields

using Gridap
import Gridap: gradient

const ConstantCellFieldLike{D,T,N,R<:FieldLike{D,T,N}} = ConstantCellValue{R}

function gradient(f::ConstantCellFieldLike)
  v = f.value
  g = gradient(v)
  ConstantCellValue(g,f.length)
end

end # module
