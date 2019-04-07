module Commons

using Numa.Helpers

export evaluate
export gradient, ∇

evaluate() = @abstractmethod
gradient() = @abstractmethod
const ∇ = gradient

end  # module Commons

# module TestModule1
# using Numa.Commons
# import Numa.Commons: evaluate, gradient
# evaluate(a::Int) = a
# end
# module TestModule2
# using Numa.Commons
# import Numa.Commons: evaluate, gradient
# evaluate(a::Real) = 2a
# end
#
# using Numa.Commons
# evaluate(1)
# evaluate(2.0)
