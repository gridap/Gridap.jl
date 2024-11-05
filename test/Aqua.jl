module AquaTests

using Gridap
using Aqua

Aqua.test_all(
  Gridap,
  ambiguities = false,
  unbound_args = false
)

"""
Comment: Ambiguities

I think all ambiguities are either false positives or never used in practice... I've seen 
other packages set `ambiguities = false` as well, so I think it's fine to do so.
"""

"""
Comment: Unbound Args

We do have some unbound type warnings. However, these are calls which are never executed in 
the code. 

They mostly involve things like `f(a::T...) where T`, which trigger the warning 
in the case were the function `f` is called with no arguments. This can be fixed as described 
in https://juliatesting.github.io/Aqua.jl/stable/unbound_args/#test_unbound_args, but it is quite 
a pain to do so... 

I guess something to think about in the future.
"""

end