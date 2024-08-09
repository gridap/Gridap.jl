"""
    ColoredArray{T,N} <: AbstractArray{T,N}

  An array with an additional color field for each entry.
    
  Note: We use a RefValue to store the array and the colors
  so that we can change the array and colors in-place without 
  making `ColoredArray` a mutable type. 
"""
struct ColoredArray{T,N} <: AbstractArray{T,N}
  array  :: Base.RefValue{Array{T,N}}
  colors :: Base.RefValue{Array{Int8,N}}
  function ColoredArray(array::Array{T,N},colors::Array{Int8,N}) where {T,N}
    new{T,N}(Base.RefValue(array),Base.RefValue(colors))
  end
end

const ColoredVector{T} = ColoredArray{T,1}
const ColoredMatrix{T} = ColoredArray{T,2}

Base.size(a::ColoredArray) = size(a.array[])

function Base.show(io::IO,o::ColoredArray)
  print(io,"ColoredArray($(o.array[]), $(o.colors[]))")
end
function Base.show(io::IO,k::MIME"text/plain",o::ColoredArray)
  print(io,"ColoredArray($(o.array[]), $(o.colors[]))")
end

function Arrays.testitem(f::ColoredArray{T}) where T
  #@notimplementedif !isconcretetype(T)
  f.array[][1], f.colors[][1]
end

"""
  ColorMap <: Map

  A map that returns a ColoredArray.

  evaluate(::ColorMap,array,colors) -> ColoredArray
"""
struct ColorMap <: Map end

function Fields.return_cache(::ColorMap,array,colors)
  return ColoredArray(array,colors)
end

function Fields.evaluate!(cache,::ColorMap,array,colors)
  cache.array[] = array
  cache.colors[] = colors
  return cache
end


# Assembly-related functions

@inline function Algebra.add_entries!(combine::Function,A::MatrixBlock,vs,is::ColoredVector,js::ColoredVector)
  _add_colored_entries!(combine,A.array,vs,is.array[],js.array[],is.colors[],js.colors[])
end

@inline function Algebra.add_entries!(combine::Function,A::VectorBlock,vs,is::ColoredVector)
  _add_colored_entries!(combine,A.array,vs,is.array[],is.colors[])
end

@inline function _add_colored_entries!(combine::Function,A,vs::Nothing,is,js,cis,cjs)
  for (lj,j) in enumerate(js)
    if j>0
      cj = cjs[lj]
      for (li,i) in enumerate(is)
        if i>0
          ci = cis[li]
          add_entry!(combine,A[ci,cj],nothing,i,j)
        end
      end
    end
  end
  A
end

@inline function _add_colored_entries!(combine::Function,A,vs,is,js,cis,cjs)
  for (lj,j) in enumerate(js)
    if j>0
      cj = cjs[lj]
      for (li,i) in enumerate(is)
        if i>0
          ci = cis[li]
          vij = vs[li,lj]
          add_entry!(combine,A[ci,cj],vij,i,j)
        end
      end
    end
  end
  A
end

@inline function _add_colored_entries!(combine::Function,A,vs::Nothing,is,cis)
  for (li, i) in enumerate(is)
    if i>0
      ci = cis[li]
      add_entry!(A[ci],nothing,i)
    end
  end
  A
end

@inline function _add_colored_entries!(combine::Function,A,vs,is,cis)
  for (li, i) in enumerate(is)
    if i>0
      ci = cis[li]
      vi = vs[li]
      add_entry!(A[ci],vi,i)
    end
  end
  A
end

struct ColorMask <: Map
  color :: Int8
end

function Fields.evaluate!(cache,k::ColorMask,v::T,c::Int8) where T <: Number
  if c == k.color
    return v
  else
    return -one(T)
  end
end
