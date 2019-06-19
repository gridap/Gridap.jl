module ConstantCellFields

using Gridap
import Gridap: gradient
import Gridap: HasGradientStyle

const ConstantCellFieldLike{D,T,N,R<:FieldLike{D,T,N}} = ConstantCellValue{R}

# TODO This returns a different object each time
function gradient(f::ConstantCellFieldLike)
  v = f.value
  @assert HasGradientStyle(f) == GradientYesStyle()
  g = gradient(v)
  ConstantCellValue(g,f.length)
end

function HasGradientStyle(::Type{ConstantCellValue{R}}) where R<:FieldLike
  HasGradientStyle(R)
end

end # module
