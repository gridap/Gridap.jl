
"""
This can be used as a CellField as long as one evaluates it
on the stored CellPoint.
"""
struct CellState{T,P<:CellPoint} <: CellField
  points::P
  values::AbstractArray

  function CellState{T}(::UndefInitializer,points::CellPoint) where T
    values = _init_values(T,get_cell_data(points))
    P = typeof(points)
    new{T,P}(points,values)
  end
  
  function CellState(v::Number,points::CellPoint)
    values = _init_values(v,get_cell_data(points))
    T = typeof(v)
    P = typeof(points)
    new{T,P}(points,values)
  end
end

function _init_values(::Type{T},x::AbstractArray{<:Point}) where T
  N = ndims(x)
  Array{T,N}(undef,size(x))
end

function _init_values(::Type{T},x::AbstractArray{<:AbstractVector{<:Point}}) where T
  [Vector{T}(undef,length(xi)) for xi in x]
end

function _init_values(v::Number,x::AbstractArray{<:Point})
  fill(v,size(x))
end

function _init_values(v::Number,x::AbstractArray{<:AbstractVector{<:Point}})
  [fill(v,length(xi)) for xi in x]
end

function get_cell_data(f::CellState)
  @unreachable """\n
  get_cell_data cannot be called on a CellState

  If you see this error messase it is likelly that you are trying to perform
  an operation that does not make sense for a CellState. In most cases,
  to do the wanted operation, you would need first to project the CellState
  to a FESpace (e.g. via a L2 projection).
  """
end

get_triangulation(f::CellState) = get_triangulation(f.points)
DomainStyle(::Type{CellState{T,P}}) where {T,P} = DomainStyle(P)

function evaluate!(cache,f::CellState,x::CellPoint)
  if f.points === x
    f.values
  else
    @unreachable """\n
    It is not possible to evaluate the given CellState on the given CellPoint.

    a CellState can only be evaluated at the CellPoint it was created from.
    If you want to evaluate at another location, you would need first to project the CellState
    to a FESpace (e.g. via a L2 projection).
    """
  end
end

function CellState{T}(::UndefInitializer,a) where T
  points = get_cell_points(a)
  CellState{T}(undef,points)
end

function CellState(v::Number,a)
  points = get_cell_points(a)
  CellState(v,points)
end

function update_state!(updater::Function,f::CellField...)
  ids = findall(map(i->isa(i,CellState),f))
  @assert length(ids) > 0 """\n
  At least one CellState object has to be given to the update_state! function
  """
  a = f[ids]
  x = first(a).points
  @assert all(map(i->i.points===x,a)) """\n
  All the CellState objects given to the update_state! function need to be
  defined on the same CellPoint.
  """
  fx = map(i->evaluate(i,x),f)
  if num_cells(x) > 0
    fxi = map(first,fx)
    fxiq = map(first,fxi)
    need_to_update, states = first_and_tail(updater(fxiq...))
    @assert isa(need_to_update,Bool) && isa(states,Tuple{Vararg{Number}}) """\n
    Wrong return value of the user-defined updater Function. The signature is

        need_to_update, state = first_and_tail(updater(args...))

    where need_to_update is a Bool telling if we need to update_state and
    states is a Tuple of Number objects with the new states.
    See `first_and_tail` for further details.
    """
    msg = """\n
    The number of new states given by the updater Function does not match
    the number of CellState objects given as arguments in the update_state! Funciton.
    """
    @check length(states) <= length(a) msg
  end
  caches = map(array_cache,fx)
  x_data = get_cell_data(x)
  cache_x = array_cache(x_data)
  _update_state_variables!(updater,caches,fx,cache_x,x_data)
  nothing
end

@noinline function  _update_state_variables!(updater,caches,fx,cache_x,x)
  ncells = length(x)
  for cell in 1:ncells
    fxi = map((c,f)->getindex!(c,f,cell),caches,fx)
    xi = getindex!(cache_x,x,cell)
    for q in 1:length(xi)
      fxiq = map(f->f[q],fxi)
      need_to_update, states = first_and_tail(updater(fxiq...))
      if need_to_update
        _update_states!(fxi,q,states,Val{length(states)}())
      end
    end
  end
end

@inline function _update_states!(b,q,states,::Val{i}) where i
  _update_state!(b,q,states,Val{i}())
  _update_states!(b,q,states,Val{i-1}())
  nothing
end

@inline function _update_states!(b,q,states,::Val{0})
  nothing
end

@inline function _update_state!(b,q,states,::Val{i}) where i
  b[i][q] = states[i]
  nothing
end

