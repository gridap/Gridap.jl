
"""
"""
macro law(fundef)
  s = "The @law macro is only allowed in function definitions"
  @assert isa(fundef,Expr) s
  @assert fundef.head in (:(=), :function) s
  funname = fundef.args[1].args[1]
  nargs = length(fundef.args[1].args)-1
  @assert nargs >= 1 "The @law macro is only allowed in functions with at leats one argument"

  y = fundef.args[1].args[2]
  x =  fundef.args[1].args[3:end]
  a = fundef.args[1].args[2:end]
  q = quote
    function $(funname)($(y)::GridapType,$(x...))
      operate($(funname),$(a...))
    end
    $(fundef)
  end
  :($(esc(q)))
end
